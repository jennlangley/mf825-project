# Load libraries
library(quadprog)
library(nloptr)
library(moments)
library(tseries)
library(ggplot2)
library(reshape2)

# Load data
data <- read.csv('./hedge_fund_data_monthly.csv')

# Fix Date and Numeric Conversion
data$Date <- as.Date(data[[1]], format = "%m/%d/%Y")

# Split data
prices <- as.matrix(sapply(data[, 2:11], as.numeric))
rf_rates <- as.numeric(data[, 12])

# Calculate returns
asset_returns <- diff(prices) / prices[-nrow(prices), ]
colnames(asset_returns) <- colnames(prices)

# Get RF returns and excess returns
rf_returns <- rf_rates[-1] / 100
excess_returns <- asset_returns - rf_returns

# Skewness and Kurtosis (after excess_returns available)
skewness_vals <- apply(excess_returns, 2, skewness)
kurtosis_vals <- apply(excess_returns, 2, kurtosis)

# Normality Check (Jarque-Bera Test)
normality_pvalues <- apply(excess_returns, 2, function(x) jarque.bera.test(x)$p.value)
normality_status <- ifelse(normality_pvalues < 0.05, "No", "Yes")

# Create summary table
normality_summary <- data.frame(
  Asset = colnames(excess_returns),
  JB_pvalue = round(normality_pvalues, 2),
  Skewness = round(skewness_vals, 2),
  Kurtosis = round(kurtosis_vals, 2),
  Normal = normality_status
)
print(normality_summary)

# Portfolio Optimization Settings
lambda <- 3
gamma <- 1
delta <- 1
gamma_power <- 3

mean_returns <- colMeans(excess_returns)
cov_matrix <- cov(excess_returns)
n_assets <- ncol(excess_returns)

# 2-Moment Portfolio (No Shorts)
Dmat <- 2 * cov_matrix
dvec <- mean_returns
Amat <- cbind(rep(1, n_assets), diag(n_assets))
bvec <- c(1, rep(0, n_assets))

opt_2 <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
weights_2 <- opt_2$solution

# 2-Moment Portfolio (With Shorts)
Amat_short <- matrix(1, n_assets, 1)
bvec_short <- 1
opt_2_short <- solve.QP(Dmat, dvec, Amat_short, bvec_short, meq = 1)
weights_2_short <- opt_2_short$solution

# Helper Functions
power_utility_fn <- function(w) {
  -(t(w) %*% mean_returns / (1 - gamma_power) - (gamma_power / 2) * t(w) %*% cov_matrix %*% w)
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

# General Optimization Function
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

# Optimize
weights_power <- optimize_portfolio(power_utility_fn, power_utility_grad)
weights_power_short <- optimize_portfolio(power_utility_fn, power_utility_grad, short_allowed = TRUE)
weights_3 <- optimize_portfolio(three_moment_fn, three_moment_grad)
weights_3_short <- optimize_portfolio(three_moment_fn, three_moment_grad, short_allowed = TRUE)
weights_4 <- optimize_portfolio(four_moment_fn, four_moment_grad)
weights_4_short <- optimize_portfolio(four_moment_fn, four_moment_grad, short_allowed = TRUE)

# Print Weights
weights_df <- data.frame(
  Asset = colnames(excess_returns),
  `2-Moment` = round(weights_2, 4),
  `Power Utility` = round(weights_power, 4),
  `3-Moment` = round(weights_3, 4),
  `4-Moment` = round(weights_4, 4),
  `2-Moment (Shorts)` = round(weights_2_short, 4),
  `Power Utility (Shorts)` = round(weights_power_short, 4),
  `3-Moment (Shorts)` = round(weights_3_short, 4),
  `4-Moment (Shorts)` = round(weights_4_short, 4)
)
print(weights_df)

# Certainty Equivalents
calculate_CE <- function(w, model = "2") {
  port_mean <- as.numeric(t(w) %*% mean_returns)
  port_var <- as.numeric(t(w) %*% cov_matrix %*% w)
  port_skew <- sum((w^3) * skewness_vals)
  port_kurt <- sum((w^4) * kurtosis_vals)
  
  if (model == "2") {
    port_mean - (lambda/2) * port_var
  } else if (model == "power") {
    port_mean / (1 - gamma_power) - (gamma_power/2) * port_var
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

# Plot Skewness and Kurtosis
sk_kurt <- data.frame(
  Asset = names(skewness_vals),
  Skewness = skewness_vals,
  Kurtosis = kurtosis_vals
)
sk_kurt_melted <- melt(sk_kurt, id.vars = "Asset")

ggplot(sk_kurt_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Asset Skewness and Kurtosis", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Portfolio Weights
weights_melted <- melt(weights_df, id.vars = "Asset")

ggplot(weights_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights Across Models", y = "Weight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


