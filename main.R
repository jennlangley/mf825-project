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

# Calculate percentage changes (returns)
asset_returns <- matrix(NA, nrow = nrow(prices) - 1, ncol = ncol(prices))
colnames(asset_returns) <- colnames(prices)

for (i in 1:ncol(prices)) {
  asset_returns[, i] <- diff(prices[, i]) / prices[1:(nrow(prices)-1), i]
}

# Get the RF returns for the same periods
rf_returns <- rf_rates[-1] / 100

# Calculate excess returns
excess_returns <- asset_returns - rf_returns

# Normality Check (Jarque-Bera Test)
normality_pvalues <- apply(excess_returns, 2, function(x) jarque.bera.test(x)$p.value)

# Report "Yes" or "No" based on a 0.05 significance level
normality_status <- ifelse(normality_pvalues < 0.05, "No", "Yes")

# Round skewness and kurtosis to 2 decimals
skewness_vals_rounded <- round(skewness_vals, 2)
kurtosis_vals_rounded <- round(kurtosis_vals, 2)

# Create a summary table with rounded values
normality_summary <- data.frame(
  Asset = colnames(excess_returns),
  JB_pvalue = round(normality_pvalues, 2),  # Rounded p-value to 2 decimals
  Skewness = skewness_vals_rounded,  # Rounded skewness
  Kurtosis = kurtosis_vals_rounded,  # Rounded kurtosis
  `Normal?` = normality_status  # "Yes" or "No" for normality
)

# Print the summary table
print(normality_summary)


# Portfolio Optimization Settings
lambda <- 3
gamma <- 1
delta <- 1
gamma_power <- 3

mean_returns <- colMeans(excess_returns)
cov_matrix <- cov(excess_returns)
skewness_vals <- apply(excess_returns, 2, skewness)
kurtosis_vals <- apply(excess_returns, 2, kurtosis)
n_assets <- ncol(excess_returns)

# 2-Moment Portfolio (Mean-Variance, no shorts)
Dmat <- 2 * cov_matrix
dvec <- mean_returns
Amat <- cbind(rep(1, n_assets), diag(n_assets))
bvec <- c(1, rep(0, n_assets))

# Without shorting
opt_2 <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
weights_2 <- opt_2$solution
cat("2-Moment Portfolio (No Shorts):\n")
print(round(weights_2, 4))

# With shorting
Amat_short <- matrix(1, nrow = n_assets, ncol = 1)
bvec_short <- 1
opt_2_short <- solve.QP(Dmat, dvec, Amat_short, bvec_short, meq = 1)
weights_2_short <- opt_2_short$solution
cat("\n2-Moment Portfolio (With Shorts):\n")
print(round(weights_2_short, 4))

# Power Utility Portfolio (no shorts)
power_utility_fn <- function(w) {
  - (t(w) %*% mean_returns / (1 - gamma_power) - (gamma_power / 2) * t(w) %*% cov_matrix %*% w)
}
power_utility_grad <- function(w) {
  - (mean_returns - gamma_power * cov_matrix %*% w)
}
opt_power <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = power_utility_fn,
  eval_grad_f = power_utility_grad,
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_power <- opt_power$solution
cat("\nPower Utility Portfolio (No Shorts):\n")
print(round(weights_power, 4))

# With shorting
opt_power_short <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = power_utility_fn,
  eval_grad_f = power_utility_grad,
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_power_short <- opt_power_short$solution
cat("\nPower Utility Portfolio (With Shorts):\n")
print(round(weights_power_short, 4))

# 3-Moment Portfolio (no shorts)
three_moment_fn <- function(w) {
  - (t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w + (gamma/6) * sum((w^3) * skewness_vals))
}
three_moment_grad <- function(w) {
  - (mean_returns - lambda * cov_matrix %*% w + (gamma/2) * (w^2) * skewness_vals)
}
opt_3 <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = three_moment_fn,
  eval_grad_f = three_moment_grad,
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_3 <- opt_3$solution
cat("\n3-Moment Portfolio (No Shorts):\n")
print(round(weights_3, 4))

# With shorting
opt_3_short <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = three_moment_fn,
  eval_grad_f = three_moment_grad,
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_3_short <- opt_3_short$solution
cat("\n3-Moment Portfolio (With Shorts):\n")
print(round(weights_3_short, 4))

# 4-Moment Portfolio (no shorts)
four_moment_fn <- function(w) {
  - (t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w +
       (gamma/6) * sum((w^3) * skewness_vals) - (delta/24) * sum((w^4) * kurtosis_vals))
}
four_moment_grad <- function(w) {
  - (mean_returns - lambda * cov_matrix %*% w + 
       (gamma/2) * (w^2) * skewness_vals - (delta/6) * (w^3) * kurtosis_vals)
}
opt_4 <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = four_moment_fn,
  eval_grad_f = four_moment_grad,
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_4 <- opt_4$solution
cat("\n4-Moment Portfolio (No Shorts):\n")
print(round(weights_4, 4))

# With shorting
opt_4_short <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = four_moment_fn,
  eval_grad_f = four_moment_grad,
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
)
weights_4_short <- opt_4_short$solution
cat("\n4-Moment Portfolio (With Shorts):\n")
print(round(weights_4_short, 4))

# Combine weights
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

# Calculate CE
calculate_CE <- function(w, mean_ret, cov_mat, skew, kurt, model = "2") {
  port_mean <- as.numeric(t(w) %*% mean_ret)
  port_var <- as.numeric(t(w) %*% cov_mat %*% w)
  port_skew <- sum((w^3) * skew)
  port_kurt <- sum((w^4) * kurt)
  
  if (model == "2") {
    return(port_mean - (lambda/2) * port_var)
  } else if (model == "power") {
    return(port_mean / (1 - gamma_power) - (gamma_power/2) * port_var)
  } else if (model == "3") {
    return(port_mean - (lambda/2) * port_var + (gamma/6) * port_skew)
  } else if (model == "4") {
    return(port_mean - (lambda/2) * port_var + (gamma/6) * port_skew - (delta/24) * port_kurt)
  }
}

CE_values <- data.frame(
  Model = c("2-Moment", "Power Utility", "3-Moment", "4-Moment"),
  CE = c(
    calculate_CE(weights_2, mean_returns, cov_matrix, skewness_vals, kurtosis_vals, "2"),
    calculate_CE(weights_power, mean_returns, cov_matrix, skewness_vals, kurtosis_vals, "power"),
    calculate_CE(weights_3, mean_returns, cov_matrix, skewness_vals, kurtosis_vals, "3"),
    calculate_CE(weights_4, mean_returns, cov_matrix, skewness_vals, kurtosis_vals, "4")
  )
)
print(CE_values)

### Plots

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

# Plot Weights Across Models
weights_melted <- melt(weights_df, id.vars = "Asset")
ggplot(weights_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights: Different Utility Models", y = "Weight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot CE Comparison
ggplot(CE_values, aes(x = Model, y = CE, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Certainty Equivalent Across Models", y = "Certainty Equivalent") +
  theme_minimal()
