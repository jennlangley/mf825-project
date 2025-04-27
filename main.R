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

# --- Utility Functions ---
# These are separated from the optimization objectives to allow separate calculation of utilities

# Expected utility calculation functions
expected_utility_2moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w
}

expected_utility_power <- function(w) {
  t(w) %*% mean_returns - (gamma_power/2) * t(w) %*% cov_matrix %*% w
}

expected_utility_3moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w + (gamma/6) * sum((w^3) * skewness_vals)
}

expected_utility_4moment <- function(w) {
  t(w) %*% mean_returns - (lambda/2) * t(w) %*% cov_matrix %*% w + 
    (gamma/6) * sum((w^3) * skewness_vals) - (delta/24) * sum((w^4) * kurtosis_vals)
}

# --- Objective Functions for Optimization ---

power_utility_fn <- function(w) {
  -expected_utility_power(w)
}
power_utility_grad <- function(w) {
  -(mean_returns - gamma_power * cov_matrix %*% w)
}

three_moment_fn <- function(w) {
  -expected_utility_3moment(w)
}
three_moment_grad <- function(w) {
  -(mean_returns - lambda * cov_matrix %*% w + (gamma/2) * (w^2) * skewness_vals)
}

four_moment_fn <- function(w) {
  -expected_utility_4moment(w)
}
four_moment_grad <- function(w) {
  -(mean_returns - lambda * cov_matrix %*% w +
      (gamma/2) * (w^2) * skewness_vals - (delta/6) * (w^3) * kurtosis_vals)
}

# --- Certainty Equivalent Calculation ---
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

# --- Optimization Function ---
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

# 2-Moment (No Shorts)
Dmat <- 2 * cov_matrix
dvec <- mean_returns
Amat <- cbind(rep(1, n_assets), diag(n_assets))
bvec <- c(1, rep(0, n_assets))
opt_2 <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
weights_2 <- opt_2$solution

# 2-Moment (With Shorts)
Amat_short <- matrix(1, n_assets, 1)
bvec_short <- 1
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

# --- Calculate Expected Utilities Across All Models and Portfolio Types ---
# For comparing utilities across models, we calculate utility using each model's own utility function

# Calculate utilities for each portfolio under its own model
utility_2_noshorting <- as.numeric(expected_utility_2moment(weights_2))
utility_power_noshorting <- as.numeric(expected_utility_power(weights_power))
utility_3_noshorting <- as.numeric(expected_utility_3moment(weights_3))
utility_4_noshorting <- as.numeric(expected_utility_4moment(weights_4))

utility_2_shorting <- as.numeric(expected_utility_2moment(weights_2_short))
utility_power_shorting <- as.numeric(expected_utility_power(weights_power_short))
utility_3_shorting <- as.numeric(expected_utility_3moment(weights_3_short))
utility_4_shorting <- as.numeric(expected_utility_4moment(weights_4_short))

# Calculate CEs for each portfolio
ce_2_noshorting <- calculate_CE(weights_2, "2")
ce_power_noshorting <- calculate_CE(weights_power, "power")
ce_3_noshorting <- calculate_CE(weights_3, "3")
ce_4_noshorting <- calculate_CE(weights_4, "4")

ce_2_shorting <- calculate_CE(weights_2_short, "2")
ce_power_shorting <- calculate_CE(weights_power_short, "power")
ce_3_shorting <- calculate_CE(weights_3_short, "3")
ce_4_shorting <- calculate_CE(weights_4_short, "4")

# Create utility and CE comparison dataframe
utility_comparison <- data.frame(
  Model = c("2-Moment", "Power Utility", "3-Moment", "4-Moment"),
  Expected_Utility_NoShorts = c(utility_2_noshorting, utility_power_noshorting, utility_3_noshorting, utility_4_noshorting),
  Expected_Utility_Shorts = c(utility_2_shorting, utility_power_shorting, utility_3_shorting, utility_4_shorting),
  CE_NoShorts = c(ce_2_noshorting, ce_power_noshorting, ce_3_noshorting, ce_4_noshorting),  
  CE_Shorts = c(ce_2_shorting, ce_power_shorting, ce_3_shorting, ce_4_shorting),
  Utility_Improvement = c(utility_2_shorting - utility_2_noshorting, 
                          utility_power_shorting - utility_power_noshorting,
                          utility_3_shorting - utility_3_noshorting,
                          utility_4_shorting - utility_4_noshorting),
  CE_Improvement = c(ce_2_shorting - ce_2_noshorting,
                     ce_power_shorting - ce_power_noshorting,
                     ce_3_shorting - ce_3_noshorting,
                     ce_4_shorting - ce_4_noshorting)
)
print(utility_comparison)

# Convert to long format for plotting
utility_long <- melt(utility_comparison[, c("Model", "Expected_Utility_NoShorts", "Expected_Utility_Shorts")], 
                     id.vars = "Model", 
                     variable.name = "Portfolio_Type", 
                     value.name = "Expected_Utility")
utility_long$Portfolio_Type <- gsub("Expected_Utility_", "", utility_long$Portfolio_Type)

ce_long <- melt(utility_comparison[, c("Model", "CE_NoShorts", "CE_Shorts")], 
                id.vars = "Model", 
                variable.name = "Portfolio_Type", 
                value.name = "Certainty_Equivalent")
ce_long$Portfolio_Type <- gsub("CE_", "", ce_long$Portfolio_Type)

# --- Plots ---

# Combined plot of Expected Utility and CE side by side
utility_plot <- ggplot(utility_long, aes(x = Model, y = Expected_Utility, fill = Portfolio_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Expected Utility by Model", y = "Expected Utility", x = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("NoShorts" = "#4DAF4A", "Shorts" = "#E41A1C"))

ce_plot <- ggplot(ce_long, aes(x = Model, y = Certainty_Equivalent, fill = Portfolio_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Certainty Equivalent by Model", y = "Certainty Equivalent", x = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("NoShorts" = "#4DAF4A", "Shorts" = "#E41A1C"))

# Improvement from allowing shorts
utility_improvement_plot <- ggplot(utility_comparison, aes(x = Model, y = Utility_Improvement)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  labs(title = "Utility Improvement from Allowing Shorts", y = "Utility Improvement", x = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ce_improvement_plot <- ggplot(utility_comparison, aes(x = Model, y = CE_Improvement)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  labs(title = "CE Improvement from Allowing Shorts", y = "CE Improvement", x = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange multiple plots
grid.arrange(utility_plot, ce_plot, utility_improvement_plot, ce_improvement_plot, ncol = 2)
ggsave("./plots/utility_ce_comparison.png", arrangeGrob(utility_plot, ce_plot, utility_improvement_plot, ce_improvement_plot, ncol = 2), width = 12, height = 10)

# Portfolio Weights Across Models (No Shorts)
weights_melted_no_shorts <- melt(weights_df[, c("Asset", "2-Moment", "Power Utility", "3-Moment", "4-Moment")], id.vars = "Asset")
colnames(weights_melted_no_shorts) <- c("Asset", "Model", "Weight")

ggplot(weights_melted_no_shorts, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights Across Models (No Shorts)", y = "Weight", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/weights_no_shorts.png", width = 8, height = 5)

# Portfolio Weights Across Models (With Shorts)
weights_melted_shorts <- melt(weights_df[, c("Asset", "2-Moment (Shorts)", "Power Utility (Shorts)", "3-Moment (Shorts)", "4-Moment (Shorts)")], id.vars = "Asset")
colnames(weights_melted_shorts) <- c("Asset", "Model", "Weight")

ggplot(weights_melted_shorts, aes(x = Asset, y = Weight, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights Across Models (With Shorts)", y = "Weight", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/weights_shorts.png", width = 8, height = 5)

# Normality Test Plot
normality_melted <- melt(normality_summary[, c("Asset", "JB_pvalue", "Shapiro_pvalue")], id.vars = "Asset")
ggplot(normality_melted, aes(x = Asset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Normality Test P-values", y = "P-value", x = "Asset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./plots/normality_tests.png", width = 10, height = 6)