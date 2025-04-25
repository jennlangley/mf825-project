# Load libraries
library(quadprog)
library(moments)
library(ggplot2)
library(tidyquant)
library(dplyr)
library(tidyr)
library(nloptr)

# Load data
tickers <- c("GLD", "SLV",     # Precious Metals
             "USO", "UNG",     # Energy
             "CORN", "WEAT",   # Agriculture
             "JO", "SOYB", "DBA")  

data <- tq_get(tickers, from = "2022-01-01", to = Sys.Date(), get = "stock.prices")

returns <- data %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted,
               mutate_fun = dailyReturn,
               col_rename = "daily_return") %>%
  pivot_wider(names_from = symbol, values_from = daily_return) %>%
  na.omit()

# Drop any columns with remaining NA values
returns <- returns %>% select(where(~ sum(is.na(.)) == 0))

returns_matrix <- as.matrix(returns %>% select(-date))

# Parameters
lambda <- 3
gamma <- 1
delta_kurtosis <- 1
gamma_power <- 3

# Stats
mean_returns <- colMeans(returns_matrix)
cov_matrix <- cov(returns_matrix)
skewness_values <- apply(returns_matrix, 2, skewness)
kurtosis_values <- apply(returns_matrix, 2, kurtosis)
n_assets <- ncol(returns_matrix)

### 1. 2-Moment Utility Portfolio (Mean-Variance Optimal)
Dmat <- 2 * cov_matrix
dvec <- mean_returns
Amat <- cbind(rep(1, n_assets), diag(n_assets))
bvec <- c(1, rep(0, n_assets))  # weights sum to 1, no short sales

opt_2 <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
weights_2 <- round(opt_2$solution, 4)

### 2. Power Utility under Lognormal
utility_power_fn <- function(w, mu, sigma, gamma) {
  port_mean <- as.numeric(t(w) %*% mu)
  port_var <- as.numeric(t(w) %*% sigma %*% w)
  - (port_mean / (1 - gamma) - (gamma / 2) * port_var)
}

opt_power <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = function(w) utility_power_fn(w, mean_returns, cov_matrix, gamma_power),
  eval_grad_f = function(w) nl.grad(w, function(w) utility_power_fn(w, mean_returns, cov_matrix, gamma_power)),
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1.0e-8)
)
weights_power <- round(opt_power$solution, 4)

### 3. 3-Moment Utility (adds Skewness)
utility_3_fn <- function(w, mu, sigma, skew, lambda, gamma) {
  - (t(w) %*% mu - (lambda / 2) * t(w) %*% sigma %*% w + (gamma / 6) * sum((w^3) * skew))
}

opt_3 <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = function(w) utility_3_fn(w, mean_returns, cov_matrix, skewness_values, lambda, gamma),
  eval_grad_f = function(w) nl.grad(w, function(w) utility_3_fn(w, mean_returns, cov_matrix, skewness_values, lambda, gamma)),
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1.0e-8)
)
weights_3 <- round(opt_3$solution, 4)

### 4. 4-Moment Utility (adds Kurtosis)
utility_4_fn <- function(w, mu, sigma, skew, kurt, lambda, gamma, delta) {
  - (t(w) %*% mu - (lambda / 2) * t(w) %*% sigma %*% w +
       (gamma / 6) * sum((w^3) * skew) -
       (delta / 24) * sum((w^4) * kurt))
}

opt_4 <- nloptr(
  x0 = rep(1/n_assets, n_assets),
  eval_f = function(w) utility_4_fn(w, mean_returns, cov_matrix, skewness_values, kurtosis_values, lambda, gamma, delta_kurtosis),
  eval_grad_f = function(w) nl.grad(w, function(w) utility_4_fn(w, mean_returns, cov_matrix, skewness_values, kurtosis_values, lambda, gamma, delta_kurtosis)),
  lb = rep(0, n_assets),
  ub = rep(1, n_assets),
  eval_g_eq = function(w) sum(w) - 1,
  eval_jac_g_eq = function(w) rep(1, n_assets),
  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1.0e-8)
)
weights_4 <- round(opt_4$solution, 4)

### Compare
comparison_df <- data.frame(
  Asset = colnames(returns_matrix),
  `2-Moment` = weights_2,
  `Power Utility` = weights_power,
  `3-Moment` = weights_3,
  `4-Moment` = weights_4
)

print(comparison_df)

### Plot Comparison
comparison_long <- pivot_longer(comparison_df, cols = -Asset, names_to = "Model", values_to = "Weight")

ggplot(comparison_long, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights Across Utility Models") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

weight_change_df <- data.frame(
  Asset = colnames(returns_matrix),
  Skewness = round(skewness_values, 3),
  `Δ3-Moment vs 2-Moment` = weights_3 - weights_2,
  `Δ4-Moment vs 2-Moment` = weights_4 - weights_2
)

# Convert to long format for plotting
weight_change_long <- pivot_longer(
  weight_change_df,
  cols = starts_with("Δ"),
  names_to = "Model_Comparison",
  values_to = "Weight_Change"
)

# Plot Skewness vs. Change in Weight
ggplot(weight_change_long, aes(x = Skewness, y = Weight_Change, label = Asset)) +
  geom_point(aes(color = Model_Comparison), size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Impact of Skewness on Portfolio Weights",
       subtitle = "Change in Allocation vs. 2-Moment Model",
       x = "Skewness of Asset",
       y = "Change in Weight",
       color = "Model Comparison") +
  theme_minimal()
