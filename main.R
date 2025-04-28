# --- Load Libraries ---
library(quadprog)
library(nloptr)
library(moments)
library(tseries)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(viridis)
library(scales)
library(dplyr) 
library(tidyr) 

# --- Load and Prepare Data ---
data <- read.csv('./data/hedge_fund_data_monthly.csv')
ff_factors <- read.csv('./data/ff_factors_monthly.csv')

# Convert the first column to Date objects
data$Date <- as.Date(data[[1]], format = "%m/%d/%Y")

# Extract asset prices (columns 2 to 11) and convert to a numeric matrix.
prices <- as.matrix(sapply(data[, 2:11], as.numeric))

# --- Risk-free Rate from Fama-French Factors ---
# Replace 'RF' below with the actual name of the risk-free rate column in ff_factors
# (e.g., "RF" or "RF..%." or similar. Check colnames(ff_factors) if unsure.)
rf_rates <- as.numeric(ff_factors$RF)  # adjust if your column is named differently

# Convert risk-free from percent to decimal and align to asset returns periods
# If RF is 0.03 (meaning 0.03%), dividing by 100 gives 0.0003 (decimal return)
rf_returns <- rf_rates[-length(rf_rates)] / 100

# --- Calculate Returns ---
# Calculate simple arithmetic asset returns.
asset_returns <- diff(prices) / prices[-nrow(prices), ]
colnames(asset_returns) <- colnames(prices) # Keep original asset names

# --- Align dates if necessary (important if the two files have different date ranges) ---
# If your files are perfectly aligned by row, you can skip this.
# Otherwise, match dates between asset_returns and rf_returns.
# Example: (uncomment and adapt if needed)
# asset_dates <- data$Date[-1]
# factor_dates <- as.Date(ff_factors$Date, format="...") # adjust date format
# idx <- match(asset_dates, factor_dates)
# rf_returns <- rf_returns[idx]

# --- Calculate excess returns (asset return - risk-free return) ---
excess_returns <- asset_returns - rf_returns

# --- Calculate Individual Asset Statistics ---
mean_returns <- colMeans(excess_returns, na.rm = TRUE)
cov_matrix <- cov(excess_returns, use = "pairwise.complete.obs")
skewness_vals <- apply(excess_returns, 2, skewness, na.rm = TRUE)

asset_stats <- data.frame(
  Asset = colnames(excess_returns),
  Mean = mean_returns,
  Variance = diag(cov_matrix),
  StdDev = sqrt(diag(cov_matrix)),
  Skewness = skewness_vals
)

print("Individual Asset Statistics (Ordered by Mean Descending):")
print(asset_stats[order(-asset_stats$Mean), ])
print("Individual Asset Statistics (Ordered by Variance Ascending):")
print(asset_stats[order(asset_stats$Variance), ])

# --- Calculate Moments and Normality Tests ---
skewness_vals <- apply(excess_returns, 2, skewness)
kurtosis_vals <- apply(excess_returns, 2, kurtosis)
normality_pvalues <- apply(excess_returns, 2, function(x) jarque.bera.test(x)$p.value)
shapiro_pvalues <- apply(excess_returns, 2, function(x) shapiro.test(x)$p.value)
normality_summary <- data.frame(
  Asset = colnames(excess_returns),
  JB_pvalue = round(normality_pvalues, 4),
  Shapiro_pvalue = round(shapiro_pvalues, 4),
  Skewness = round(skewness_vals, 3),
  Kurtosis = round(kurtosis_vals, 3),
  JB_Normal = ifelse(normality_pvalues < 0.05, "No", "Yes"),
  Shapiro_Normal = ifelse(shapiro_pvalues < 0.05, "No", "Yes")
)
print("Normality Summary:")
print(normality_summary)

# --- Optimization Parameters ---
lambda <- 3 # Risk aversion for variance
gamma <- 1  # Skewness preference
delta <- 1  # Kurtosis aversion (kept for potential future use)
mean_returns <- colMeans(excess_returns)
cov_matrix <- cov(excess_returns)
n_assets <- ncol(excess_returns)
n_obs <- nrow(excess_returns)

# --- Coskewness Tensor ---
coskewness_tensor <- function(R) {
  N <- ncol(R)
  T <- nrow(R)
  mu <- colMeans(R)
  S <- array(0, dim = c(N, N, N))
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  for (i in 1:N) {
    for (j in 1:N) {
      for (k in 1:N) {
        S[i, j, k] <- mean(R_centered[, i] * R_centered[, j] * R_centered[, k])
      }
    }
  }
  S
}
S_tensor <- coskewness_tensor(excess_returns)

portfolio_variance <- function(w, Sigma) {
  as.numeric(t(w) %*% Sigma %*% w)
}
portfolio_coskewness <- function(w, S) {
  N <- length(w)
  skew_val <- 0
  for(i in 1:N) {
    for(j in 1:N) {
      for(k in 1:N) {
        skew_val <- skew_val + w[i]*w[j]*w[k]*S[i,j,k]
      }
    }
  }
  return(skew_val)
}

expected_utility_2moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * portfolio_variance(w, cov_matrix)
}
expected_utility_3moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * portfolio_variance(w, cov_matrix) +
    (gamma/6) * portfolio_coskewness(w, S_tensor)
}

three_moment_grad <- function(w) {
  N <- length(w)
  grad_coskew_term <- numeric(N)
  for (i in 1:N) {
    sum_jk = 0
    for (j in 1:N) {
      for (k in 1:N) {
        sum_jk = sum_jk + w[j] * w[k] * S_tensor[i, j, k]
      }
    }
    grad_coskew_term[i] = sum_jk
  }
  gradient <- -(mean_returns - lambda * cov_matrix %*% w + (gamma/2) * grad_coskew_term)
  return(gradient)
}
three_moment_fn <- function(w) { -expected_utility_3moment(w) }

optimize_portfolio_nloptr <- function(fn, grad) {
  nloptr(
    x0 = rep(1/n_assets, n_assets),
    eval_f = fn,
    eval_grad_f = grad,
    lb = rep(0, n_assets),
    ub = rep(1, n_assets),
    eval_g_eq = function(w) sum(w) - 1,
    eval_jac_g_eq = function(w) rep(1, n_assets),
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8, maxeval = 1000)
  )$solution
}

# --- Perform Optimizations (Long-Only) ---
Dmat_2m <- 2 * lambda * cov_matrix
Amat_2m <- cbind(rep(1, n_assets), diag(n_assets))
bvec_2m <- c(1, rep(0, n_assets))
opt_2 <- solve.QP(Dmat = lambda * cov_matrix, dvec = mean_returns, Amat = Amat_2m, bvec = bvec_2m, meq = 1)
weights_2 <- opt_2$solution
weights_2[weights_2 < 1e-6] <- 0
weights_2 <- weights_2 / sum(weights_2)

weights_3 <- optimize_portfolio_nloptr(three_moment_fn, three_moment_grad)
weights_3[weights_3 < 1e-6] <- 0
weights_3 <- weights_3 / sum(weights_3)

weights_df_longonly <- data.frame(
  Asset = colnames(excess_returns),
  `2-Moment (MV)` = round(weights_2, 4),
  `3-Moment (MVS)` = round(weights_3, 4),
  check.names = FALSE
)
print("Optimal Long-Only Weights:")
print(weights_df_longonly)

portfolio_stats <- function(w) {
  port_mean <- as.numeric(t(w) %*% mean_returns)
  port_var <- portfolio_variance(w, cov_matrix)
  port_skew <- portfolio_coskewness(w, S_tensor)
  return(data.frame(Mean = port_mean, Variance = port_var, Skewness = port_skew))
}
stats_2 <- portfolio_stats(weights_2)
stats_3 <- portfolio_stats(weights_3)
stats_compare <- rbind(stats_2, stats_3)
stats_compare$Model <- c("2-Moment (MV)", "3-Moment (MVS)")
print("Portfolio Characteristics Comparison:")
print(stats_compare)

ce_2 <- expected_utility_2moment(weights_2)
ce_3 <- expected_utility_3moment(weights_3)
ce_df <- data.frame(
  Model = c("2-Moment (MV)", "3-Moment (MVS)"),
  CertaintyEquivalent = c(ce_2, ce_3)
)
print("Certainty Equivalent Comparison:")
print(ce_df)

# --- Plots ---
weights_melted <- melt(weights_df_longonly, id.vars = "Asset", variable.name = "Model", value.name = "Weight")
plot_weights <- ggplot(weights_melted, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Optimal Long-Only Portfolio Weights",
       subtitle = "Comparing Mean-Variance (MV) vs. Mean-Variance-Skewness (MVS)",
       y = "Weight", x = "Asset") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis_d(option = "D")
print(plot_weights)
ggsave("./plots/weights_longonly_comparison.png", plot_weights, width = 10, height = 6)

stats_melted <- stats_compare %>%
  select(Model, Mean, Variance, Skewness) %>%
  pivot_longer(cols = c(Mean, Variance, Skewness), names_to = "Statistic", values_to = "Value")
plot_stats <- ggplot(stats_melted, aes(x = Statistic, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~Statistic, scales = "free_y") +
  labs(title = "Portfolio Characteristics Comparison",
       subtitle = "MV-Optimal vs. MVS-Optimal Portfolios (Long-Only)",
       y = "Value", x = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis_d(option = "D")
print(plot_stats)
ggsave("./plots/portfolio_stats_comparison.png", plot_stats, width = 10, height = 5)

plot_ce <- ggplot(ce_df, aes(x = Model, y = CertaintyEquivalent, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.5f", CertaintyEquivalent)), vjust = -0.5) +
  labs(title = "Certainty Equivalent Comparison (Long-Only)",
       subtitle = "Calculated using each model's own utility function",
       y = "Certainty Equivalent", x = "Model") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_viridis_d(option = "D")
print(plot_ce)
ggsave("./plots/ce_longonly_comparison.png", plot_ce, width = 7, height = 5)

normality_melted <- melt(normality_summary[, c("Asset", "JB_pvalue", "Shapiro_pvalue")], id.vars = "Asset")
plot_normality <- ggplot(normality_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = 0.05, label = "p = 0.05", hjust = 1.1, vjust = -0.5, color = "red") +
  labs(title = "Normality Test P-values for Asset Excess Returns", y = "P-value", x = "Asset") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "top") +
  scale_fill_brewer(palette = "Set1", labels = c("Jarque-Bera", "Shapiro-Wilk"))
print(plot_normality)
ggsave("./plots/normality_tests.png", plot_normality, width = 10, height = 6)