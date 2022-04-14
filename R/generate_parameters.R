#' Pick parameter values for simulation, from a normal distribution
#'
#' @param Tt Final period (t=0,1,...,Tt)
#' @param mu_Beta_W1 numeric. mean for parameters in covariate models
#' @param mu_Beta_W2 numeric. mean for parameters in covariate models
#' @param mu_Beta_A numeric. mean for parameters in treatment models
#' @param mu_Beta_Y numeric. mean for parameters in outcome models
#' @param sd_Beta_W1 numeric. standard deviation for parameters in covariate models
#' @param sd_Beta_W2 numeric. standard deviation for parameters in covariate models
#' @param sd_Beta_A numeric. standard deviation for parameters in treatment models
#' @param sd_Beta_Y numeric. standard deviation for parameters in outcome models
#' @param Beta_YU numeric. constant coefficient for effect of unmeasured, time-invariant, binary covariate U in Y models
#' @param range_ymeans numeric vector of length 2. minimum and maximum values for marginal means of Y
#'
#' @return list of matrices, one for each variable to be generated (U, W1, W2, A, Y). Each matrix has T+1 rows and a different number of columns corresponding to the coefficients of a (generalized) linear model. The columns names refer to the term in that model.
#' @export
#'
#'
#' @examples
#' generate_parameters(Tt=2)
#'
#' @family simulation functions
generate_parameters <- function(Tt,
                                mu_Beta_W1=0.2, mu_Beta_W2=0.2, mu_Beta_A=0.2, mu_Beta_Y=0.2,
                                sd_Beta_W1=0.2, sd_Beta_W2=0.2, sd_Beta_A=0.2, sd_Beta_Y=0.2,
                                range_ymeans = c(-5, 5),
                                Beta_YU = 0.1) {

  Beta_U = -log(1 / 0.5 - 1)
  mean_U = stats::plogis(Beta_U)

  # matrices of parameters
  Beta_W1 = matrix( stats::rnorm(n = 2*(Tt+1), mean =mu_Beta_W1, sd = sd_Beta_W1), nrow = Tt + 1, byrow = TRUE)# 2 because terms for A + intercept
  Beta_W2 = matrix( stats::rnorm(n = 2*(Tt+1), mean =mu_Beta_W2, sd = sd_Beta_W2), nrow = Tt + 1, byrow = TRUE)# 2 because terms for A + intercept
  Beta_A = matrix( stats::rnorm(n = 5*(Tt+1), mean =mu_Beta_A, sd = sd_Beta_A), nrow = Tt + 1, byrow = TRUE)# 5 because terms for U, L, W, W^2, + intercept
  Beta_Y = matrix( stats::rnorm(n = 6*(Tt+1), mean =mu_Beta_Y, sd = sd_Beta_Y), nrow = Tt + 1, byrow = TRUE)# 6 because terms for U, L, W, W^2, A + intercept

  Beta_Y[,2] = Beta_YU #constant effect of U0 = part of sufficient conditions for parallel trends

  # balancing intercepts
  mean_Y = stats::runif(Tt + 1, min = range_ymeans[1], max = range_ymeans[2])
  mean_W1 = stats::runif(Tt + 1, min = 0.1, max = 0.5)
  mean_W2 = 0
  mean_A_conditional = stats::runif(Tt + 1, min = 1/(Tt^2), max = 1/(Tt^2) + 1/Tt) # this is a function of T to avoid very small pr(A_t=0) for large t
  mean_A_marginal = 1 - cumprod(1 - mean_A_conditional)  # basically a kaplan meier risk function

  Beta_Y[, 1] = mean_Y - Beta_Y[, 2]*mean_U - Beta_Y[, 3]*mean_W1 - (Beta_Y[,4] + Beta_Y[,5])*mean_W2 - Beta_Y[,6]*mean_A_marginal
  Beta_A[, 1] = -log(1 / mean_A_conditional - 1) - Beta_A[, 2]*mean_U - Beta_A[, 3]*mean_W1 - (Beta_A[, 4] + Beta_A[,5])*mean_W2
  Beta_W1[, 1] = -log(1 / mean_W1 - 1) - Beta_W1[, 2]*mean_A_marginal #fix - should depend on At-1 not At
  Beta_W2[, 1] = mean_W2 - Beta_W2[, 2]*mean_A_marginal #fix - should depend on At-1 not At

  colnames(Beta_Y) = c('intercept', 'U0', 'W1t', 'W2t','W2t^2','At')
  colnames(Beta_A) = c('intercept', 'U0', 'W1t', 'W2t','Wt^2')
  colnames(Beta_W1) = colnames(Beta_W2) = c('intercept', 'At-1')

  return ( list(U=Beta_U, W1=Beta_W1, W2=Beta_W2, A=Beta_A, Y=Beta_Y))

}
