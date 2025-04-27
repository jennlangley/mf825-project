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

# --- Load and Prepare Data ---
data <- read.csv('./hedge_fund_data_monthly.csv')
data$Date <- as.Date(data[[1]], format = "%m/%d/%Y")
prices <- as.matrix(sapply(data[, 2:11], as.numeric))
rf_rates <- as.numeric(data[, 12])

asset_returns <- diff(prices) / prices[-nrow(prices), ]
colnames(asset_returns) <- colnames(prices)
rf_returns <- rf_rates[-1] / 100
excess_returns <- asset_returns - rf_returns

# --- Calculate Moments and Normality Tests ---
skewness_vals <- apply(excess_returns, 2, skewness)
kurtosis_vals <- apply(excess_returns, 2, kurtosis)

# Jarque-Bera Test
normality_pvalues <- apply(excess_returns, 2, function(x) jarque.bera.test(x)$p.value)
normality_status <- ifelse(normality_pvalues < 0.05, "No", "Yes")

# Shapiro-Wilk Test
shapiro_pvalues <- apply(excess_returns, 2, function(x) shapiro.test(x)$p.value)
shapiro_status <- ifelse(shapiro_pvalues < 0.05, "No", "Yes")

# Normality Summary Table
normality_summary <- data.frame(
  Asset = colnames(excess_returns),
  JB_pvalue = round(normality_pvalues, 4),
  Shapiro_pvalue = round(shapiro_pvalues, 4),
  Skewness = round(skewness_vals, 3),
  Kurtosis = round(kurtosis_vals, 3),
  JB_Normal = normality_status,
  Shapiro_Normal = shapiro_status
)
print(normality_summary)

# --- Optimization Parameters ---
lambda <- 3
gamma <- 1
delta <- 1
gamma_power <- 3

mean_returns <- colMeans(excess_returns)
cov_matrix <- cov(excess_returns)
n_assets <- ncol(excess_returns)

# --- Objective Functions ---

# Corrected Power Utility Function
power_utility_fn <- function(w) {
  -(t(w) %*% mean_returns - (gamma_power / 2) * t(w) %*% cov_matrix %*% w)
}
power_utility_grad <- function(w) {
  -(mean_returns - gamma_power * cov_matrix %*% w)
}

three_moment_fn <- function(w) {
  -(t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w + (gamma/6) * sum((w^3) * skewness_vals))
}
three_moment_grad <- function(w) {
  -(mean_returns - lambda * cov_matrix %*% w + (gamma/2) * (w^2) * skewness_vals)
}

four_moment_fn <- function(w) {
  -(t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w +
      (gamma/6) * sum((w^3) * skewness_vals) - (delta/24) * sum((w^4) * kurtosis_vals))
}
four_moment_grad <- function(w) {
  -(mean_returns - lambda * cov_matrix %*% w +
      (gamma/2) * (w^2) * skewness_vals - (delta/6) * (w^3) * kurtosis_vals)
}

# --- General Optimization Function ---
optimize_portfolio <- function(fn, grad, short_allowed = FALSE) {
  nloptr(
    x0 = rep(1/n_assets, n_assets),
    eval_f = fn,
    eval_grad_f = grad,
    lb = if (short_allowed) rep(-1, n_assets) else rep(0, n_assets),
    ub = rep(1, n_assets),
    eval_g_eq = function(w) sum(w) - 1,
    eval_jac_g_eq = function(w) rep(1, n_assets),
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
  )$solution
}

# --- Solve for Weights ---
# 2-Moment (No Shorts)
Dmat <- 2 * cov_matrix
dvec <- mean_returns
Amat <- cbind(rep(1, n_assets), diag(n_assets))
bvec <- c(1, rep(0, n_assets))
opt_2 <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
weights_2 <- opt_2$solution

# 2-Moment (With Shorts)
Amat_short <- matrix(1, n_assets, 1)
opt_2_short <- solve.QP(Dmat, dvec, Amat_short, bvec_short, meq = 1)
weights_2_short <- opt_2_short$solution

# Other models
weights_power <- optimize_portfolio(power_utility_fn, power_utility_grad)
weights_power_short <- optimize_portfolio(power_utility_fn, power_utility_grad, short_allowed = TRUE)
weights_3 <- optimize_portfolio(three_moment_fn, three_moment_grad)
weights_3_short <- optimize_portfolio(three_moment_fn, three_moment_grad, short_allowed = TRUE)
weights_4 <- optimize_portfolio(four_moment_fn, four_moment_grad)
weights_4_short <- optimize_portfolio(four_moment_fn, four_moment_grad, short_allowed = TRUE)

# --- Portfolio Weights Table ---
weights_df <- data.frame(
  Asset = colnames(excess_returns),
  `2-Moment` = round(weights_2, 4),
  `Power Utility` = round(weights_power, 4),
  `3-Moment` = round(weights_3, 4),
  `4-Moment` = round(weights_4, 4),
  `2-Moment (Shorts)` = round(weights_2_short, 4),
  `Power Utility (Shorts)` = round(weights_power_short, 4),
  `3-Moment (Shorts)` = round(weights_3_short, 4),
  `4-Moment (Shorts)` = round(weights_4_short, 4),
  check.names = FALSE
)
print(weights_df)

# --- Certainty Equivalents ---
calculate_CE <- function(w, model = "2") {
  port_mean <- as.numeric(t(w) %*% mean_returns)
  port_var <- as.numeric(t(w) %*% cov_matrix %*% w)
  port_skew <- sum((w^3) * skewness_vals)
  port_kurt <- sum((w^4) * kurtosis_vals)
  
  if (model == "2") {
    port_mean - (lambda/2) * port_var
  } else if (model == "power") {
    port_mean - (gamma_power/2) * port_var
  } else if (model == "3") {
    port_mean - (lambda/2) * port_var + (gamma/6) * port_skew
  } else if (model == "4") {
    port_mean - (lambda/2) * port_var + (gamma/6) * port_skew - (delta/24) * port_kurt
  }
}

CE_values <- data.frame(
  Model = c("2-Moment", "Power Utility", "3-Moment", "4-Moment"),
  CE = c(
    calculate_CE(weights_2, "2"),
    calculate_CE(weights_power, "power"),
    calculate_CE(weights_3, "3"),
    calculate_CE(weights_4, "4")
  )
)
print(CE_values)

# --- Plots ---

# Portfolio Weights Across Models (No Shorts)
weights_melted_no_shorts <- melt(weights_df[, c("Asset", "2-Moment", "Power Utility", "3-Moment", "4-Moment")], id.vars = "Asset")
colnames(weights_melted_no_shorts) <- c("Asset", "Model", "Weight")

ggplot(weights_melted_no_shorts, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights Across Models (No Shorts)", y = "Weight", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/weights_no_shorts.png", width = 8, height = 5)

# Certainty Equivalents
ggplot(CE_values, aes(x = Model, y = CE, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Certainty Equivalents Comparison", y = "Certainty Equivalent", x = "Model") +
  theme_minimal()
ggsave("./plots/certainty_equivalents.png", width = 6, height = 4)

# Change in Weights vs 2-Moment
weights_diff <- data.frame(Asset = weights_df$Asset)
weights_diff$`Power Utility Diff` <- weights_df$`Power Utility` - weights_df$`2-Moment`
weights_diff$`3-Moment Diff` <- weights_df$`3-Moment` - weights_df$`2-Moment`
weights_diff$`4-Moment Diff` <- weights_df$`4-Moment` - weights_df$`2-Moment`
weights_diff_melted <- melt(weights_diff, id.vars = "Asset")
colnames(weights_diff_melted) <- c("Asset", "Model", "Weight_Difference")

ggplot(weights_diff_melted, aes(x = Asset, y = Weight_Difference, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Change in Portfolio Weights vs 2-Moment", y = "Weight Difference", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/change_in_weights.png", width = 8, height = 5)

# Normality Test Plot
normality_melted <- melt(normality_summary[, c("Asset", "JB_pvalue", "Shapiro_pvalue")], id.vars = "Asset")
ggplot(normality_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Normality Test P-values", y = "P-value", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/normality_tests.png", width = 10, height = 6)
