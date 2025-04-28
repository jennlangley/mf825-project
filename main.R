# Portfolio Optimization with Higher Order Moments
library(nloptr)
library(moments)

# Import and process data
data <- read.csv('./data/hedge_fund_data_monthly.csv')
ff_factors <- read.csv('./data/ff_factors_monthly.csv')
data$Date <- as.Date(data[[1]], format = "%m/%d/%Y")
prices <- as.matrix(sapply(data[, 2:11], as.numeric))
rf_returns <- as.numeric(ff_factors$RF[-1]) / 100

# Calculate returns
asset_returns <- diff(log(prices))
colnames(asset_returns) <- colnames(prices)

# Remove Fund 3; bad
asset_returns <- asset_returns[, -3]  
excess_returns <- asset_returns - rf_returns

# Annualized statistics
mean_returns <- colMeans(asset_returns) * 12
variance_vals <- (round(100 * apply(asset_returns, 2, sd) * sqrt(12), 2))^2
cov_matrix <- cov(asset_returns) * 12
skewness_vals <- apply(asset_returns, 2, skewness)
kurtosis_vals <- apply(asset_returns, 2, kurtosis)
standard_deviation <- sqrt(diag(cov_matrix))
sharpe_ratios <- mean_returns / standard_deviation
rho <- apply(asset_returns, 2, function(x) {
  cor(x[-1], x[-length(x)])
})

asset_stats <- data.frame(
  Asset     = colnames(asset_returns),
  Mean      = round(colMeans(asset_returns)*12*100, 2),
  Variance  = variance_vals,
  StdDev    = round(apply(asset_returns,2,sd)*sqrt(12)*100, 2),
  Sharpe    = round(mean_returns/standard_deviation, 2),
  Skewness  = round(skewness_vals, 2),
  Kurtosis  = round(kurtosis_vals, 2),
  Autocorrelation = round(rho, 3)
)

# --- Normality Test ---
# Compute Shapiro–Wilk
shapiro_p <- apply(asset_returns, 2, function(x) shapiro.test(x)$p.value)

asset_stats$Shapiro_p <- round(shapiro_p, 4)
asset_stats$Normal   <- ifelse(shapiro_p < 0.05, "No", "Yes")


print("Individual Asset Statistics:") # ordered by Mean Desc
print(asset_stats[order(-asset_stats$Mean), ])

# QQ plots
oldpar <- par(mfrow = c(3, 3),            #3x3 grid
              mar   = c(2, 2, 2, 1))     # margin

for(i in seq_len(ncol(asset_returns))) {
  # extract the series
  series <- asset_returns[, i]
  
  # compute line parameters from 1st and 3rd quartiles
  qy <- quantile(series, c(0.25, 0.75))
  zx <- qnorm(c(0.25, 0.75))
  slope     <- diff(qy) / diff(zx)
  intercept <- qy[1] - slope * zx[1]
    qq <- qqnorm(series, plot.it = FALSE)
    plot(qq,
       main   = colnames(asset_returns)[i],
       xlab   = "", ylab = "",
       xlim   = range(qq$x), 
       ylim   = range(qq$y),
       xaxt   = "n",  # no x-axis labels
       yaxt   = "n",  # no y-axis labels
       frame.plot = TRUE)
    abline(intercept, slope, col = "red", lwd = 1.5)
}

# restore plotting parameters
par(oldpar)


# Coskewness and cokurtosis tensors
coskewness_tensor <- function(R) {
  N <- ncol(R); S <- array(0, dim = c(N, N, N))
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  for (i in 1:N) for (j in 1:N) for (k in 1:N) 
    S[i, j, k] <- mean(R_centered[, i] * R_centered[, j] * R_centered[, k])
  S
}
cokurtosis_tensor <- function(R) {
  N <- ncol(R); K <- array(0, dim = c(N, N, N, N))
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  for (i in 1:N) for (j in 1:N) for (k in 1:N) for (l in 1:N) 
    K[i, j, k, l] <- mean(R_centered[, i] * R_centered[, j] * R_centered[, k] * R_centered[, l])
  K
}
S_tensor <- coskewness_tensor(excess_returns)
K_tensor <- cokurtosis_tensor(excess_returns)

# Portfolio moment functions
portfolio_variance <- function(w, Sigma) as.numeric(t(w) %*% Sigma %*% w)
portfolio_coskewness <- function(w, S) sum(outer(outer(w, w), w) * S)
portfolio_cokurtosis <- function(w, K) sum(outer(outer(outer(w, w), w), w) * K)

# Utility setup
lambda <- 1.0 # risk aversion
gamma <- 2.0 # skewness preference
delta <- 0.5 # kurtosis aversion

expected_utility <- list(
  "2-Moment" = function(w) sum(w * mean_returns) - (lambda/2) * portfolio_variance(w, cov_matrix),
  "3-Moment" = function(w) sum(w * mean_returns) - (lambda/2) * portfolio_variance(w, cov_matrix) + (gamma/6) * portfolio_coskewness(w, S_tensor),
  "4-Moment" = function(w) sum(w * mean_returns) - (lambda/2) * portfolio_variance(w, cov_matrix) + (gamma/6) * portfolio_coskewness(w, S_tensor) - (delta/24) * portfolio_cokurtosis(w, K_tensor)
)

# Optimization function
optimize_portfolio <- function(util_fn) {
  N <- ncol(excess_returns)
  eval_f <- function(w) {
    return(-util_fn(w))
  }
  eval_g_eq <- function(w) {
    return(sum(w) - 1)
  }
  res <- nloptr::nloptr(
    x0 = rep(1/N, N),
    eval_f = eval_f,
    lb = rep(0, N),
    ub = rep(1, N),
    eval_g_eq = eval_g_eq,
    opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-8)
  )
  res$solution
}

# Optimize all three portfolios
optimized_weights <- lapply(expected_utility, optimize_portfolio)

# Portfolio Weights Table
weights_df <- do.call(cbind, lapply(optimized_weights, function(w) round(w, 4)))
weights_df <- data.frame(Asset = colnames(excess_returns), weights_df)
print("Optimal Portfolio Weights:")
print(weights_df)

# Portfolio stats
portfolio_stats <- function(w) {
  mean_ret <- sum(w * mean_returns)
  var_ret <- portfolio_variance(w, cov_matrix)
  data.frame(
    Mean = mean_ret * 100,
    StdDev = sqrt(var_ret) * 100,
    Sharpe = mean_ret / sqrt(var_ret),
    Skewness = portfolio_coskewness(w, S_tensor),
    Kurtosis = portfolio_cokurtosis(w, K_tensor)
  )
}

stats_compare <- do.call(rbind, lapply(optimized_weights, portfolio_stats))
stats_compare$Model <- names(expected_utility)
print("Portfolio Characteristics:")
print(stats_compare)

# Utility and Certainty Equivalent Table
utility_table <- sapply(optimized_weights, function(w) sapply(expected_utility, function(fn) round(fn(w), 6)))
CE <- diag(utility_table)  # Certainty equivalent under own model

utility_comparison <- data.frame(
  Portfolio = names(expected_utility),
  MV_Utility = utility_table[1, ],
  MVS_Utility = utility_table[2, ],
  MVSK_Utility = utility_table[3, ],
  CertaintyEquivalent = CE
)
print("Expected Utility and Certainty Equivalent:")
print(utility_comparison)

model_colors <- c("steelblue", "forestgreen", "firebrick")
model_names <- c("2-Moment", "3-Moment", "4-Moment")

# Convert to matrix for barplot
weights_matrix <- t(as.matrix(weights_df[, -1]))
colnames(weights_matrix) <- weights_df$Asset

# grouped barplot
barplot(
  weights_matrix,
  beside = TRUE,
  col = model_colors,
  ylim = c(0, 1),
  ylab = "Portfolio Weight",
  main = "Optimal Portfolio Allocations by Moment Model",
  las = 2,
  cex.names = 0.9
)
legend("topright", legend = model_names, fill = model_colors, bty = "n")

# Sharpe Ratio
bp1 <- barplot(
  stats_compare$Sharpe,
  names.arg = stats_compare$Model,
  col = model_colors,
  ylab = "Sharpe Ratio",
  main = "Sharpe Ratio by Model"
)
text(
  x = bp1,
  y = stats_compare$Sharpe,
  labels = round(stats_compare$Sharpe, 2),
  pos = 3, cex = 1, offset = 0.2
)

# Skewness
bp2 <- barplot(
  stats_compare$Skewness,
  names.arg = stats_compare$Model,
  col = model_colors,
  ylab = "Skewness",
  main = "Portfolio Skewness"
)
text(
  x = bp2,
  y = stats_compare$Skewness,
  labels = signif(stats_compare$Skewness, 4),
  pos = 3, cex = 1, offset = 0.2
)
# Kurtosis
bp3 <- barplot(
  stats_compare$Kurtosis,
  names.arg = stats_compare$Model,
  col = model_colors,
  ylab = "Kurtosis",
  main = "Portfolio Kurtosis"
)
text(
  x = bp3,
  y = stats_compare$Kurtosis,
  labels = signif(stats_compare$Kurtosis, 2),
  pos = 3, cex = 1, offset = 0.2
)
# Certainty Equivalent
bp4 <- barplot(
  utility_comparison$CertaintyEquivalent,
  names.arg = utility_comparison$Portfolio,
  col = model_colors,
  ylab = "Certainty Equivalent",
  main = "Certainty Equivalent by Model"
)
text(
  x = bp4,
  y = utility_comparison$CertaintyEquivalent,
  labels = signif(utility_comparison$CertaintyEquivalent, 4),
  pos = 3, cex = 1, offset = 0.2
)