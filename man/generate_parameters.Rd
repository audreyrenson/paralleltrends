% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_parameters.R
\name{generate_parameters}
\alias{generate_parameters}
\title{Pick parameter values for simulation, from a normal distribution}
\usage{
generate_parameters(
  Tt,
  mu_Beta_W1 = 0.2,
  mu_Beta_W2 = 0.2,
  mu_Beta_A = 0.2,
  mu_Beta_Y = 0.2,
  sd_Beta_W1 = 0.2,
  sd_Beta_W2 = 0.2,
  sd_Beta_A = 0.2,
  sd_Beta_Y = 0.2,
  range_ymeans = c(-5, 5),
  Beta_YU = 0.1
)
}
\arguments{
\item{Tt}{Final period (t=0,1,...,Tt)}

\item{mu_Beta_W1}{numeric. mean for parameters in covariate models}

\item{mu_Beta_W2}{numeric. mean for parameters in covariate models}

\item{mu_Beta_A}{numeric. mean for parameters in treatment models}

\item{mu_Beta_Y}{numeric. mean for parameters in outcome models}

\item{sd_Beta_W1}{numeric. standard deviation for parameters in covariate models}

\item{sd_Beta_W2}{numeric. standard deviation for parameters in covariate models}

\item{sd_Beta_A}{numeric. standard deviation for parameters in treatment models}

\item{sd_Beta_Y}{numeric. standard deviation for parameters in outcome models}

\item{range_ymeans}{numeric vector of length 2. minimum and maximum values for marginal means of Y}

\item{Beta_YU}{numeric. constant coefficient for effect of unmeasured, time-invariant, binary covariate U in Y models}
}
\value{
list of matrices, one for each variable to be generated (U, W1, W2, A, Y). Each matrix has T+1 rows and a different number of columns corresponding to the coefficients of a (generalized) linear model. The columns names refer to the term in that model.
}
\description{
Pick parameter values for simulation, from a normal distribution
}
\examples{
generate_parameters(Tt=2)

}
\seealso{
Other simulation functions: 
\code{\link{generate_data}()}
}
\concept{simulation functions}
