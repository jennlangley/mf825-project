# ==================================================
# PORTFOLIO OPTIMIZATION WITH HIGHER-ORDER MOMENTS
# ==================================================

# --- Load essential libraries ---
library(quadprog)    # For quadratic programming
library(nloptr)      # For nonlinear optimization
library(moments)     # For statistical moments
library(tseries)     # For time series analysis
library(ggplot2)     # For plotting
library(dplyr)       # For data manipulation
library(scales)      # For pretty scales in plots

# --- Import and prepare data ---
data <- read.csv('./data/hedge_fund_data_monthly.csv')
ff_factors <- read.csv('./data/ff_factors_monthly.csv')

# Convert dates and extract numeric data
data$Date <- as.Date(data[[1]], format = "%m/%d/%Y")
prices <- as.matrix(sapply(data[, 2:11], as.numeric))
rf_rates <- as.numeric(ff_factors$RF) / 100  # Risk-free rates
rf_returns <- rf_rates[-length(rf_rates)]

# Calculate returns and excess returns
asset_returns <- diff(prices) / prices[-nrow(prices), ]
colnames(asset_returns) <- colnames(prices)
excess_returns <- asset_returns - rf_returns

# --- Calculate basic statistics for each asset ---
mean_returns <- colMeans(excess_returns, na.rm = TRUE)
cov_matrix <- cov(excess_returns, use = "pairwise.complete.obs")
skewness_vals <- apply(excess_returns, 2, skewness, na.rm = TRUE)
kurtosis_vals <- apply(excess_returns, 2, kurtosis, na.rm = TRUE)

# Calculate Sharpe ratios for individual assets
sharpe_ratios <- mean_returns / sqrt(diag(cov_matrix))

# Compile asset statistics in a data frame
asset_stats <- data.frame(
  Asset = colnames(excess_returns),
  Mean = mean_returns,
  Variance = diag(cov_matrix),
  StdDev = sqrt(diag(cov_matrix)),
  Sharpe = sharpe_ratios,
  Skewness = skewness_vals,
  Kurtosis = kurtosis_vals
)

# Display statistics sorted by Mean return
print("Individual Asset Statistics (Ordered by Mean Descending):")
print(asset_stats[order(-asset_stats$Mean), ])

# --- Test for normality (Shapiro-Wilk only) ---
shapiro_pvalues <- apply(excess_returns, 2, function(x) shapiro.test(x)$p.value)

# Compile normality test results
normality_summary <- data.frame(
  Asset = colnames(excess_returns),
  Shapiro_pvalue = round(shapiro_pvalues, 4),
  Skewness = round(skewness_vals, 3),
  Kurtosis = round(kurtosis_vals, 3),
  Normal = ifelse(shapiro_pvalues < 0.05, "No", "Yes")
)
print("Normality Summary (Shapiro-Wilk Test):")
print(normality_summary)

# --- Set optimization parameters ---
# Risk aversion, skewness preference, and kurtosis aversion parameters
lambda <- 1.0   # Risk aversion for variance (moderate)
gamma <- 2.0   # Skewness preference (strong)
delta <- 0.5    # Kurtosis aversion (low)
n_assets <- ncol(excess_returns)
n_obs <- nrow(excess_returns)

# --- Higher-order moment tensors ---
# Calculate coskewness tensor
coskewness_tensor <- function(R) {
  N <- ncol(R)
  S <- array(0, dim = c(N, N, N))
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  
  for (i in 1:N) {
    for (j in 1:N) {
      for (k in 1:N) {
        S[i, j, k] <- mean(R_centered[, i] * R_centered[, j] * R_centered[, k], na.rm = TRUE)
      }
    }
  }
  return(S)
}

# Calculate cokurtosis tensor
cokurtosis_tensor <- function(R) {
  N <- ncol(R)
  K <- array(0, dim = c(N, N, N, N))
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  
  for (i in 1:N) {
    for (j in 1:N) {
      for (k in 1:N) {
        for (l in 1:N) {
          K[i, j, k, l] <- mean(R_centered[, i] * R_centered[, j] * 
                                  R_centered[, k] * R_centered[, l], na.rm = TRUE)
        }
      }
    }
  }
  return(K)
}

# Calculate the tensors (this may take some time)
S_tensor <- coskewness_tensor(excess_returns)
K_tensor <- cokurtosis_tensor(excess_returns)

# --- Portfolio moment functions ---
# Functions to calculate portfolio moments given weights
portfolio_variance <- function(w, Sigma) {
  as.numeric(t(w) %*% Sigma %*% w)
}

portfolio_coskewness <- function(w, S) {
  N <- length(w)
  val <- 0
  for(i in 1:N) {
    for(j in 1:N) {
      for(k in 1:N) {
        val <- val + w[i]*w[j]*w[k]*S[i,j,k]
      }
    }
  }
  return(val)
}

portfolio_cokurtosis <- function(w, K) {
  N <- length(w)
  val <- 0
  for(i in 1:N) {
    for(j in 1:N) {
      for(k in 1:N) {
        for(l in 1:N) {
          val <- val + w[i]*w[j]*w[k]*w[l]*K[i,j,k,l]
        }
      }
    }
  }
  return(val)
}

# Function to calculate Sharpe ratio
portfolio_sharpe <- function(w) {
  return(sum(w * mean_returns) / sqrt(portfolio_variance(w, cov_matrix)))
}

# Function to compile portfolio statistics
portfolio_stats <- function(weights) {
  mean_ret <- sum(weights * mean_returns)
  var_ret <- portfolio_variance(weights, cov_matrix)
  skew_ret <- portfolio_coskewness(weights, S_tensor)
  kurt_ret <- portfolio_cokurtosis(weights, K_tensor)
  sharpe_ret <- mean_ret / sqrt(var_ret)
  
  data.frame(
    Mean = mean_ret,
    Variance = var_ret,
    StdDev = sqrt(var_ret),
    Sharpe = sharpe_ret,
    Skewness = skew_ret,
    Kurtosis = kurt_ret
  )
}

# --- Utility functions for different models ---
# Two-moment utility (mean-variance)
expected_utility_2moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * portfolio_variance(w, cov_matrix)
}

# Three-moment utility (mean-variance-skewness)
expected_utility_3moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * portfolio_variance(w, cov_matrix) +
    (gamma/6) * portfolio_coskewness(w, S_tensor)
}

# Four-moment utility (mean-variance-skewness-kurtosis)
expected_utility_4moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * portfolio_variance(w, cov_matrix) +
    (gamma/6) * portfolio_coskewness(w, S_tensor) - (delta/24) * portfolio_cokurtosis(w, K_tensor)
}

# --- Gradient functions for optimization ---
# Gradient for three-moment utility
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

# Function to minimize for three-moment model
three_moment_fn <- function(w) { -expected_utility_3moment(w) }

# Gradient for four-moment utility
four_moment_grad <- function(w) {
  N <- length(w)
  grad_coskew_term <- numeric(N)
  grad_cokurt_term <- numeric(N)
  
  # Coskewness part
  for (i in 1:N) {
    sum_jk = 0
    for (j in 1:N) for (k in 1:N)
      sum_jk = sum_jk + w[j] * w[k] * S_tensor[i, j, k]
    grad_coskew_term[i] = sum_jk
  }
  
  # Cokurtosis part
  for (i in 1:N) {
    sum_jkl = 0
    for (j in 1:N) for (k in 1:N) for (l in 1:N)
      sum_jkl = sum_jkl + w[j] * w[k] * w[l] * K_tensor[i, j, k, l]
    grad_cokurt_term[i] = sum_jkl
  }
  
  gradient <- -(mean_returns - lambda * cov_matrix %*% w +
                  (gamma/2) * grad_coskew_term - (delta/6) * grad_cokurt_term)
  return(gradient)
}

# Function to minimize for four-moment model
four_moment_fn <- function(w) { -expected_utility_4moment(w) }

# --- Portfolio optimization function ---
# Unconstrained optimization (will normalize weights later)
optimize_portfolio_unconstrained <- function(fn, grad) {
  result <- nloptr(
    x0 = rep(1/n_assets, n_assets),  # Equal weight initial guess
    eval_f = fn,                      # Objective function
    eval_grad_f = grad,               # Gradient function
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8, maxeval = 1000)
  )
  
  # Return normalized weights
  return(result$solution)
}

# --- Perform optimizations ---
# Two-moment optimization (mean-variance)
weights_2_unconstrained <- optimize_portfolio_unconstrained(
  function(w) -expected_utility_2moment(w),
  function(w) -(mean_returns - lambda * cov_matrix %*% w)
)
weights_2_unconstrained <- weights_2_unconstrained / sum(weights_2_unconstrained)

# Three-moment optimization (mean-variance-skewness)
weights_3_unconstrained <- optimize_portfolio_unconstrained(three_moment_fn, three_moment_grad)
weights_3_unconstrained <- weights_3_unconstrained / sum(weights_3_unconstrained)

# Four-moment optimization (mean-variance-skewness-kurtosis)
weights_4_unconstrained <- optimize_portfolio_unconstrained(four_moment_fn, four_moment_grad)
weights_4_unconstrained <- weights_4_unconstrained / sum(weights_4_unconstrained)

# --- Create summary tables ---
# Portfolio weights table
weights_df_unconstrained <- data.frame(
  Asset = colnames(excess_returns),
  `2-Moment` = round(weights_2_unconstrained, 4),
  `3-Moment` = round(weights_3_unconstrained, 4),
  `4-Moment` = round(weights_4_unconstrained, 4),
  check.names = FALSE
)
print("Optimal Unconstrained Weights:")
print(weights_df_unconstrained)

# Calculate portfolio characteristics
stats_2_unconstrained <- portfolio_stats(weights_2_unconstrained)
stats_3_unconstrained <- portfolio_stats(weights_3_unconstrained)
stats_4_unconstrained <- portfolio_stats(weights_4_unconstrained)

stats_compare_unconstrained <- rbind(stats_2_unconstrained, stats_3_unconstrained, stats_4_unconstrained)
stats_compare_unconstrained$Model <- c("2-Moment", "3-Moment", "4-Moment")
print("Portfolio Characteristics Comparison (Unconstrained):")
print(stats_compare_unconstrained)

# --- Calculate expected utility and certainty equivalents ---
# For 2-moment (MV) portfolio
mv_mean <- sum(weights_2_unconstrained * mean_returns)
mv_variance <- portfolio_variance(weights_2_unconstrained, cov_matrix)
mv_skewness <- portfolio_coskewness(weights_2_unconstrained, S_tensor)
mv_kurtosis <- portfolio_cokurtosis(weights_2_unconstrained, K_tensor)
mv_sharpe <- mv_mean / sqrt(mv_variance)

# Calculate utility under different models
mv_utility_2m <- expected_utility_2moment(weights_2_unconstrained)
mv_utility_3m <- expected_utility_3moment(weights_2_unconstrained)
mv_utility_4m <- expected_utility_4moment(weights_2_unconstrained)

# For 3-moment portfolio
mvs_mean <- sum(weights_3_unconstrained * mean_returns)
mvs_variance <- portfolio_variance(weights_3_unconstrained, cov_matrix)
mvs_skewness <- portfolio_coskewness(weights_3_unconstrained, S_tensor)
mvs_kurtosis <- portfolio_cokurtosis(weights_3_unconstrained, K_tensor)
mvs_sharpe <- mvs_mean / sqrt(mvs_variance)

# Calculate utility under different models
mvs_utility_2m <- expected_utility_2moment(weights_3_unconstrained)
mvs_utility_3m <- expected_utility_3moment(weights_3_unconstrained)
mvs_utility_4m <- expected_utility_4moment(weights_3_unconstrained)

# For 4-moment portfolio
mvsk_mean <- sum(weights_4_unconstrained * mean_returns)
mvsk_variance <- portfolio_variance(weights_4_unconstrained, cov_matrix)
mvsk_skewness <- portfolio_coskewness(weights_4_unconstrained, S_tensor)
mvsk_kurtosis <- portfolio_cokurtosis(weights_4_unconstrained, K_tensor)
mvsk_sharpe <- mvsk_mean / sqrt(mvsk_variance)

# Calculate utility under different models
mvsk_utility_2m <- expected_utility_2moment(weights_4_unconstrained)
mvsk_utility_3m <- expected_utility_3moment(weights_4_unconstrained)
mvsk_utility_4m <- expected_utility_4moment(weights_4_unconstrained)

# Compile comprehensive utility comparison
utility_comparison <- data.frame(
  Portfolio = c("2-Moment", "3-Moment", "4-Moment"),
  Mean = c(mv_mean, mvs_mean, mvsk_mean),
  Variance = c(mv_variance, mvs_variance, mvsk_variance),
  Sharpe = c(mv_sharpe, mvs_sharpe, mvsk_sharpe),
  `2M_Utility` = c(mv_utility_2m, mvs_utility_2m, mvsk_utility_2m),
  `3M_Utility` = c(mv_utility_3m, mvs_utility_3m, mvsk_utility_3m),
  `4M_Utility` = c(mv_utility_4m, mvs_utility_4m, mvsk_utility_4m)
)

print("Comprehensive Utility Comparison:")
print(utility_comparison)

# Display a more focused expected utility table
ce_df_unconstrained <- data.frame(
  Model = c("2-Moment", "3-Moment", "4-Moment"),
  Sharpe = round(c(mv_sharpe, mvs_sharpe, mvsk_sharpe), 4),
  MV_Utility = round(c(mv_utility_2m, mvs_utility_2m, mvsk_utility_2m), 6),
  MVS_Utility = round(c(mv_utility_3m, mvs_utility_3m, mvsk_utility_3m), 6),
  MVSK_Utility = round(c(mv_utility_4m, mvs_utility_4m, mvsk_utility_4m), 6),
  CertaintyEquivalent = round(c(mv_utility_2m, mvs_utility_3m, mvsk_utility_4m), 6),
  check.names = FALSE
)
print("Expected Utility and Certainty Equivalent Comparison:")
print(ce_df_unconstrained)

# --- Create visualizations ---
# Set a consistent theme for all plots
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_blank()
  )

# 1. Portfolio weights comparison
weights_long <- weights_df_unconstrained %>%
  tidyr::pivot_longer(cols = -Asset, names_to = "Model", values_to = "Weight")

plot_weights <- ggplot(weights_long, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Optimal Portfolio Weights",
    subtitle = "Comparison across different utility models",
    y = "Weight", x = ""
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set1") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_weights)
ggsave("./plots/portfolio_weights_comparison.png", plot_weights, width = 10, height = 6)

# 2. Expected utility comparison across all models
ce_long <- ce_df_unconstrained %>%
  tidyr::pivot_longer(cols = c(MV_Utility, MVS_Utility, MVSK_Utility), 
                      names_to = "UtilityType", values_to = "Utility")

plot_utility <- ggplot(ce_long, aes(x = Model, y = Utility, fill = UtilityType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Expected Utility Comparison",
    subtitle = "Portfolio performance under different utility specifications",
    y = "Expected Utility", x = ""
  ) +
  scale_fill_brewer(palette = "Set2", 
                    labels = c("Mean-Variance", "Mean-Variance-Skewness", "Mean-Variance-Skewness-Kurtosis")) +
  my_theme

print(plot_utility)
ggsave("./plots/expected_utility_comparison_all.png", plot_utility, width = 10, height = 6)

# 3. Certainty equivalent comparison
plot_ce <- ggplot(ce_df_unconstrained, aes(x = reorder(Model, CertaintyEquivalent), 
                                           y = CertaintyEquivalent, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.6f", CertaintyEquivalent)), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold") +
  labs(
    title = "Certainty Equivalent Comparison",
    subtitle = "Higher values indicate better performance given investor preferences",
    y = "Certainty Equivalent", x = ""
  ) +
  scale_fill_brewer(palette = "Set1") +
  my_theme +
  theme(legend.position = "none") +
  coord_flip()

print(plot_ce)
ggsave("./plots/certainty_equivalent_comparison.png", plot_ce, width = 9, height = 5)

# 4. Sharpe ratio comparison
plot_sharpe <- ggplot(ce_df_unconstrained, aes(x = reorder(Model, Sharpe), 
                                               y = Sharpe, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.4f", Sharpe)), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold") +
  labs(
    title = "Sharpe Ratio Comparison",
    subtitle = "Risk-adjusted return across different portfolio models",
    y = "Sharpe Ratio", x = ""
  ) +
  scale_fill_brewer(palette = "Set1") +
  my_theme +
  theme(legend.position = "none") +
  coord_flip()

print(plot_sharpe)
ggsave("./plots/sharpe_ratio_comparison.png", plot_sharpe, width = 9, height = 5)

# 5. Portfolio characteristics comparison
stats_long <- stats_compare_unconstrained %>%
  tidyr::pivot_longer(cols = c(Mean, Variance, Sharpe, Skewness, Kurtosis), 
                      names_to = "Statistic", values_to = "Value")

plot_stats <- ggplot(stats_long, aes(x = Statistic, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Statistic, scales = "free_y") +
  labs(
    title = "Portfolio Characteristics",
    subtitle = "Comparison across different utility models",
    y = "Value", x = ""
  ) +
  scale_fill_brewer(palette = "Set1") +
  my_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

print(plot_stats)
ggsave("./plots/portfolio_characteristics.png", plot_stats, width = 10, height = 6)