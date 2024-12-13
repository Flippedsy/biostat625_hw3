% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M_gaussian.R
\name{M_gaussian}
\alias{M_gaussian}
\title{Generalized Linear Model for Gaussian Family}
\usage{
M_gaussian(formula, data)
}
\arguments{
\item{formula}{A symbolic description of the model to be fitted.}

\item{data}{A data frame containing the variables in the model.}
}
\value{
A list containing coefficients, residual standard error, R-squared, adjusted R-squared,
        F-statistic, p-value, null deviance, residual deviance, and AIC.
}
\description{
This function performs a Gaussian generalized linear model using ordinary least squares (OLS).
}
\examples{
library(haven)
data <- data("completedata", package = "GLM.g.p")
result <- M_gaussian(Depression ~ Fatalism, data = data)

}
