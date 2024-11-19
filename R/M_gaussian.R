#' Generalized Linear Model for Gaussian Family
#'
#' This function performs a Gaussian generalized linear model using ordinary least squares (OLS).
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#'
#' @return A list containing coefficients, residual standard error, R-squared, adjusted R-squared,
#'         F-statistic, p-value, null deviance, residual deviance, and AIC.
#'
#' @examples
#' library(haven)
#' data <- read_sas(system.file("data", "completedata.sas7bdat", package = "GLM.g.p"))
#' result <- M_gaussian(Depression ~ Fatalism, data = data)
#'
#' @export
M_gaussian <- function(formula, data) {
  # Parse formula
  terms <- all.vars(formula)
  response <- terms[1]
  predictors <- terms[-1]

  # Extract response and predictors
  Y <- data[[response]]
  X <- as.matrix(cbind(1, data[predictors]))  # Add intercept column
  colnames(X) <- c("(Intercept)", predictors)

  # Compute coefficients using Ordinary Least Squares
  beta <- solve(t(X) %*% X, t(X) %*% Y)

  # Predicted values and residuals
  Yhat <- X %*% beta
  residuals <- Y - Yhat

  # Residual summary
  residual_summary <- quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1))

  # Residual sum of squares and mean squared error
  SSE <- sum(residuals^2)
  MSE <- SSE / (nrow(X) - ncol(X))

  # Residual standard error
  residual_standard_error <- sqrt(MSE)

  # Coefficients standard error and t-values
  XTX_inv <- solve(t(X) %*% X)
  se_beta <- sqrt(diag(XTX_inv * MSE))
  t_values <- beta / se_beta

  # P-values
  p_values <- 2 * (1 - pt(abs(t_values), df = nrow(X) - ncol(X)))

  # Null deviance
  null_Yhat <- mean(Y)
  null_deviance <- sum((Y - null_Yhat)^2)

  # Residual deviance
  residual_deviance <- SSE

  # R-squared and adjusted R-squared
  r_squared <- 1 - SSE / null_deviance
  adj_r_squared <- 1 - (1 - r_squared) * (nrow(X) - 1) / (nrow(X) - ncol(X))

  # F-statistic
  MSR <- (null_deviance - SSE) / (ncol(X) - 1)
  F_statistic <- MSR / MSE
  f_p_value <- 1 - pf(F_statistic, df1 = ncol(X) - 1, df2 = nrow(X) - ncol(X))

  # AIC calculation
  log_likelihood <- -0.5 * (nrow(X) * log(2 * pi * MSE) + SSE / MSE)
  AIC <- -2 * log_likelihood + 2 * ncol(X)

  # Output results
  cat("Deviance Residuals: \n")
  cat(sprintf("%10s %10s %10s %10s %10s\n", "Min", "1Q", "Median", "3Q", "Max"))
  cat(sprintf("%10.4f %10.4f %10.4f %10.4f %10.4f\n\n", residual_summary[1], residual_summary[2], residual_summary[3], residual_summary[4], residual_summary[5]))

  coef_table <- data.frame(
    Estimate = beta,
    `Std. Error` = se_beta,
    `t value` = t_values,
    `Pr(>|t|)` = p_values
  )
  rownames(coef_table) <- colnames(X)

  cat("Coefficients:\n")
  printCoefmat(coef_table, signif.stars = TRUE)

  cat("\n(Dispersion parameter for gaussian family taken to be", round(MSE, 6), ")\n\n")

  cat("    Null deviance:", round(null_deviance, 2), " on", nrow(X) - 1, " degrees of freedom\n")
  cat("Residual deviance:", round(residual_deviance, 2), " on", nrow(X) - ncol(X), " degrees of freedom\n")
  cat("AIC:", round(AIC, 2), "\n\n")

  cat("Number of Fisher Scoring iterations: 1\n")

  # Return results as a list (optional)
  list(
    coefficients = coef_table,
    residual_standard_error = residual_standard_error,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    f_statistic = F_statistic,
    f_p_value = f_p_value,
    null_deviance = null_deviance,
    residual_deviance = residual_deviance,
    AIC = AIC
  )
}
