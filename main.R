# Load libraries
library(quadprog)
library(moments)
library(ggplot2)
library(tidyquant)
library(dplyr)
library(tidyr)  

# Load data
tickers <- c("GLD", "SLV",    # Precious Metals
             "USO", "UNG",    # Energy
             "CORN", "WEAT",  # Agriculture
             "JO", "SOYB", "COW")  # Soft commodity and Agriculture

data <- tq_get(tickers, from = "2022-01-01", to = Sys.Date(), get = "stock.prices")

available_symbols <- unique(data$symbol)
cat("Successfully downloaded data for:", paste(available_symbols, collapse = ", "), "\n")

returns <- data %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted,
               mutate_fun = dailyReturn,
               col_rename = "daily_return") %>%
  pivot_wider(names_from = symbol, values_from = daily_return) %>%
  na.omit()  # removes dates with missing values

cat("Dimensions of returns data:", dim(returns), "\n")

# Check skewness
skewness <- returns %>%
  select(-date) %>%
  summarise(across(everything(), skewness))
print(skewness)

returns_matrix <- as.matrix(returns %>% select(-date))  # Exclude date column