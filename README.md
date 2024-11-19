# GLM.g.p: Generalized Linear Models in R

The `GLM.g.p` package provides custom implementations of Generalized Linear Models (GLMs) for Gaussian and Poisson distributions. These functions (`M_gaussian` and `M_poisson`) are designed to perform linear regression and Poisson regression using standard statistical methods such as Ordinary Least Squares (OLS) and Iteratively Reweighted Least Squares (IRLS).

## Features

- **`M_gaussian`**: Perform linear regression for Gaussian-distributed response variables.
- **`M_poisson`**: Fit Poisson regression models for count data using IRLS.

## Installation

To install the development version of `GLM.g.p` from GitHub, use the following commands in R:

```{r}
# Install devtools if not already installed
install.packages("devtools")

# Install GLM.g.p from GitHub
devtools::install_github("yourusername/GLM.g.p")
```

Usage
After installing the package, load it using:

```{r}
library(GLM.g.p)
```
Example 1: Linear Regression with M_gaussian
The M_gaussian function fits a linear model to Gaussian-distributed data. Below is an example:

```{r}
# Example data
library(haven)
data <- read_sas(system.file("data", "completedata.sas7bdat", package = "GLM.g.p"))

# Fit a Gaussian regression model
result_gaussian <- M_gaussian(Depression ~ Fatalism, data = data)

# View results
result_gaussian
```
Example 2: Poisson Regression with M_poisson
The M_poisson function fits a Poisson regression model for count data. Below is an example:

```{r}
# Example data
data <- read_sas(system.file("data", "completedata.sas7bdat", package = "GLM.g.p"))

# Fit a Poisson regression model
result_poisson <- M_poisson(Depression ~ Comorbidity1, data = data)

# View results
cat("Poisson Regression Results:\n")
cat("Coefficients:\n")
print(result_poisson$coefficients)
cat("\nNull Deviance:", result_poisson$null_deviance, "\n")
cat("Residual Deviance:", result_poisson$residual_deviance, "\n")
cat("Converged:", result_poisson$converged, "\n")
```

Functions
M_gaussian
Description: Performs linear regression using Ordinary Least Squares (OLS) for Gaussian-distributed data.
Usage:

```{r}
M_gaussian(formula, data)
```
Arguments:
formula: A symbolic description of the model to be fitted.
data: A data frame containing the variables in the model.
M_poisson
Description: Fits a Poisson regression model for count data using Iteratively Reweighted Least Squares (IRLS).
Usage:

```{r}
M_poisson(formula, data)
```

Arguments:
formula: A symbolic description of the model to be fitted.
data: A data frame containing the variables in the model.
Vignettes
For detailed examples and explanations, see the vignette:

```{r}
vignette("glm-usage", package = "GLM.g.p")
```
License
This package is licensed under the MIT License. See LICENSE for detail

