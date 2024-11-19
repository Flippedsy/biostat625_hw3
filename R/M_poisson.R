#' Poisson Regression with IRLS
#'
#' This function performs Poisson regression using Iteratively Reweighted Least Squares (IRLS).
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#'
#' @return A list containing coefficients, null deviance, residual deviance,
#'         the number of iterations, and convergence status.
#'
#' @examples
#' library(haven)
#' data <- read_sas("completedata.sas7bdat")
#' model_poisson <- M_poisson( Comorbidity1~ Depression, data = data)
#' cat("\nPoisson Regression:\n")
#' cat("Coefficients:\n")
#' print(model_poisson$coefficients)
#' cat("Null Deviance:", model_poisson$null_deviance, "\n")
#' cat("Residual Deviance:", model_poisson$residual_deviance, "\n")
#' cat("Converged:", model_poisson$converged, "\n")
#'
#' @export
M_poisson <- function(formula, data) {
  # Parse formula
  terms <- all.vars(formula)
  response <- terms[1]
  predictors <- terms[-1]

  Y <- data[[response]]
  X <- as.matrix(cbind(1, data[predictors]))  # Add intercept column

  # Initialize variables
  mu <- rep(mean(Y), length(Y))
  link_function <- function(mu) log(mu)  # Log link
  inverse_link_function <- function(eta) exp(eta)

  beta <- rep(0, ncol(X))
  eta <- link_function(mu)

  # Iteratively Reweighted Least Squares (IRLS)
  max_iter <- 100
  tol <- 1e-6
  converged <- FALSE

  for (iter in 1:max_iter) {
    weights <- mu
    W <- diag(as.vector(weights))
    z <- eta + (Y - mu) / mu

    beta_new <- solve(t(X) %*% W %*% X, t(X) %*% W %*% z)
    eta <- X %*% beta_new
    mu <- inverse_link_function(eta)

    if (sum(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }

  if (!converged) {
    warning("IRLS did not converge.")
  }

  # Deviance calculations
  log_likelihood <- sum(Y * log(mu) - mu - lfactorial(Y))
  null_mu <- mean(Y)
  null_log_likelihood <- sum(Y * log(null_mu) - null_mu - lfactorial(Y))

  null_deviance <- -2 * null_log_likelihood
  residual_deviance <- -2 * log_likelihood

  list(
    coefficients = beta,
    null_deviance = null_deviance,
    residual_deviance = residual_deviance,
    iterations = iter,
    converged = converged
  )
}
