#' Append n lags of a group of variables to a dataset
#'
#' @param data data frame
#' @param n_lags int. number of lags to append
#' @param lag_vars character vector. Names of dataframe columns for which you want lags
#' @param default vector of length 1 or nrow(data), as in the `default` argument in `dplyr::lag`.
#'
#' @return The original data frame with additional columns `x_lag1`,`x_lag2`,...,`x_lag{n_lags}` for `x` in `lag_vars`.
#' @export
#'
#' @examples
#' expand.grid(id=1:5, tt=1:3) %>%
#'    dplyr::mutate(x=rnorm(15)) %>%
#'    dplyr::group_by(id) %>%
#'    append_lags(n_lags = 2, lag_vars='x')
append_lags = function(data, n_lags, lag_vars, default=NA) {
  for(n in 1:n_lags)
    for(lag_var in lag_vars)
      #this uses mutate because you can't take advantage of dplyr's groups otherwise
      data[[glue::glue('{lag_var}_lag{n}')]] = dplyr::mutate(data, v=dplyr::lag(.data[[lag_var]], n, default))$v
  return (data)
}



safe_predict = function(...) {
  withCallingHandlers({
    stats::predict(...)
  }, warning = function(w) {
    if (startsWith(conditionMessage(w), "prediction from a rank-deficient fit may be misleading"))
      invokeRestart("muffleWarning")
  })
}


make_time_dummies = function (df, timevar)
{
  new_timevars = paste0(timevar, min(df[[timevar]]):max(df[[timevar]]))
  time_columns = stats::model.matrix(stats::as.formula(glue::glue("~-1 +factor({timevar})")),
                              data = df)
  colnames(time_columns) = new_timevars
  result = dplyr::bind_cols(df, tibble::as_tibble(time_columns))
  attr(result, "timevars") = new_timevars
  return(result)
}


dens <- function(x, model, newdata, suppress_rank_warnings=FALSE, binomial_n=1) {

  #Evaluate data against the implied density from a fitted `glm` object based on possibly new covariate values, with automatic detection of family and link.
  #
  # @param x vector of values at which to evaluate the density
  # @param model a fitted `glm` object
  # @param newdata data frame with values of the covariates based on which to simulate the response, similar to `predict.glm`. nrow(newdata) must equal length(x).
  # @param binomial_n int length nrow(newdata). If `model` was fit with `link=binomial`, you can optionally pass a vector of group sizes to simulate aggregate binomial data.

  if(length(x) != nrow(newdata)) stop('length(x) must equal nrow(newdata)')
  if(!length(binomial_n) %in% c(1, length(x))) stop('binomial n must be length 1 or length(x)')

  pred_fun = if(suppress_rank_warnings) safe_predict else stats::predict

  N = length(x)
  n = binomial_n
  eta = pred_fun(model, newdata=newdata, type='link')
  dispersion = summary(model)$dispersion
  dens_fun = get_densfun(stats::family(model))

  dens_fun(x, eta, dispersion, n)
}

get_densfun <- function(family) {
  # Generate a function that looks up a density in accordance with the specification in a `stats::family` object
  #
  # @description This function is useful for evaluating data against the implied density of a fitted `glm` object
  #
  # @param family an object of type `stats::family`
  #
  # @return A function with arguments `x` (data), `eta`, (linear predictor from `glm` object), `dispersion` (`summary(m)$dispersion` where `m` is a fitted `glm`), and  `n`, vector of group sizes (optional unless family=binomial)
  # @export
  #
  # @examples
  #
  # fam = Gamma('log')
  # densfun = get_densfun(fam)
  # densfun(0.2, eta = 2, dispersion=0.5)
  #
  # # the above is equivalent to:
  # dgamma(0.2, scale = exp(2), shape = 1/0.5)
  function(x, eta, dispersion, n) {
    mu = family$linkinv(eta)
    stopifnot(length(x) == length(eta))

    switch ( family$family,
             'binomial' = stats::dbinom(x, size=n, prob=mu),
             'gaussian' = stats::dnorm(x, mu, sqrt(dispersion)),
             'Gamma' = stats::dgamma(x, shape = 1/dispersion, scale = dispersion*mu),
             'poisson' = stats::dpois(x, lambda = mu))
  }
}
