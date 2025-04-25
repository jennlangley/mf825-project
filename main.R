# Load libraries
library(quadprog)
library(moments)
library(ggplot2)
library(tidyquant)
library(dplyr)
library(tidyr)  

# Load data
tickers <- c("GLD", "SLV", "USO", "UNG", "CORN", "WEAT", "JO", "SOYB", "COW")
data <- tq_get(tickers, from = "2022-01-01", to = Sys.Date(), get = "stock.prices")

available_symbols <- unique(data$symbol)
cat("Successfully downloaded data for:", paste(available_symbols, collapse = ", "), "\n")

returns <- data %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted,
               mutate_fun = dailyReturn,
               col_rename = "daily_return") %>%
  pivot_wider(names_from = symbol, values_from = daily_return) %>%
  na.omit()  # Removes dates with missing values

# Convert returns to matrix for calculations
returns_matrix <- as.matrix(returns %>% select(-date))  # Exclude date column

# Calculate mean, variance, skewness, kurtosis
mean_returns <- colMeans(returns_matrix)
variance_returns <- apply(returns_matrix, 2, var)
skewness_values <- apply(returns_matrix, 2, skewness)
kurtosis_values <- apply(returns_matrix, 2, kurtosis)

# Define risk aversion and skewness preference parameters
lambda <- 3  # Risk aversion
gamma <- 1   # Skewness preference
delta_kurtosis <- 1  # Kurtosis preference
gamma_power <- 3  # Risk aversion for power utility (lognormal)

# Define a function to calculate utility
calculate_utility <- function(mean, variance, skewness = NULL, kurtosis = NULL, utility_type = "2-moment") {
  if (utility_type == "2-moment") {
    return(mean - (lambda / 2) * variance)
  } else if (utility_type == "3-moment") {
    return(mean - (lambda / 2) * variance + (gamma / 6) * skewness)
  } else if (utility_type == "power") {
    return(mean / (1 - gamma_power) - (gamma_power / 2) * variance)
  } else if (utility_type == "4-moment") {
    return(mean - (lambda / 2) * variance + (gamma / 6) * skewness - (delta_kurtosis / 24) * kurtosis)
  }
}

# Calculate utilities for all types
utility_2_moment <- calculate_utility(mean_returns, variance_returns, utility_type = "2-moment")
utility_3_moment <- calculate_utility(mean_returns, variance_returns, skewness_values, utility_type = "3-moment")
utility_power <- calculate_utility(mean_returns, variance_returns, utility_type = "power")
utility_4_moment <- calculate_utility(mean_returns, variance_returns, skewness_values, kurtosis_values, utility_type = "4-moment")

# Data frame for comparison
utilities_comparison_all <- data.frame(
  Symbol = colnames(returns_matrix),
  Utility_2_Moment = utility_2_moment,
  Utility_3_Moment = utility_3_moment,
  Utility_Power = utility_power,
  Utility_4_Moment = utility_4_moment
)

# Print comparison
print(utilities_comparison_all)

# Plot comparison
utilities_comparison_all_long <- utilities_comparison_all %>%
  pivot_longer(cols = starts_with("Utility"), names_to = "Utility_Type", values_to = "Utility_Value")

ggplot(utilities_comparison_all_long, aes(x = Symbol, y = Utility_Value, fill = Utility_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Different Utility Functions",
       x = "Asset",
       y = "Utility Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
