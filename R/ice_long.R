ice_pipeline_long <- function(df_obs,
                              df_interv,
                              inside_formula_t,
                              inside_formula_tmin1,
                              outside_formula,
                              Tt,
                              t_col,
                              n_nested=Tt,
                              inside_family='gaussian',
                              pt_link_fun=NULL,
                              binomial_n=NULL,
                              tibble=TRUE,
                              models=TRUE,
                              tmle=FALSE,
                              weights=NULL,
                              suppress_rank_warnings=FALSE
                              ) {

  ice_preds = recursive_ice_long(Tt=Tt,
                                 n_nested=n_nested,
                                 df_obs = df_obs,
                                 df_interv=df_interv,
                                 inside_formula_t = inside_formula_t,
                                 inside_formula_tmin1 = inside_formula_tmin1,
                                 outside_formula = outside_formula,
                                 inside_family=inside_family,
                                 t_col = t_col,
                                 binomial_n = binomial_n,
                                 models=models,
                                 tmle=tmle,
                                 weights=weights,
                                 suppress_rank_warnings=suppress_rank_warnings)

  ice_estimates = estimate_ice_long(ice_preds, t_col, pt_link_fun, binomial_n)[-1] #[-1] because 1st entry will be NA for tmin1, by design

  if(!tibble) {
    result = ice_estimates
  } else {
    result = tibble::tibble(t=1:Tt,
                            estimate = ice_estimates)
  }

  if(models) attr(result, 'models') = attr(ice_preds, 'models')

  return(result)

}

recursive_ice_long <- function(Tt,
                               n_nested,
                               df_obs,
                               df_interv,
                               inside_formula_t,
                               inside_formula_tmin1,
                               outside_formula,
                               inside_family,
                               t_col,
                               binomial_n=NULL,
                               models=TRUE,
                               tmle=FALSE,
                               weights=NULL,
                               suppress_rank_warnings=FALSE) {

  ### Adding convenience variables, so that users can reference
  ### {tvars} and {n} in formulas.
  tvars = attr(df_obs, 'timevars')
  tvars = paste0('(', paste( tvars[max(1, n_nested) : length(tvars)] , collapse='+'), ')') #this is like (t{n} + t{n+1} + ... + t{Tt}) for referencing formulas
  n = n_nested #for convenience, allow user to specify formulas with x{n} where {n} is shorthand for {n_nested}

  if(n_nested == 0) {
    ## innermost Q models
    two_models = fit_inner_models_long(df_obs=df_obs,
                                        formula_t = glue::glue(inside_formula_t),
                                        formula_tmin1 = glue::glue(inside_formula_tmin1),
                                        family=inside_family)

    preds_n = sapply(X = two_models,
                     FUN = if(suppress_rank_warnings) safe_predict else predict,
                     newdata = df_interv,
                     type = 'response',
                     simplify = TRUE)

    if(tmle) {
      if(is.null(weights)) stop("Must specify weights if tmle==TRUE")
      for(s in c('tmin1','t')) {
        preds_n[, s] = tmle_update(preds_n[, s],
                                   y = df_obs$Y,
                                   w = df_obs[[glue::glue(weights)]],
                                   family=gaussian,
                                   subset=df_obs$A==0)
      }
    }

    if(models) attr(preds_n, 'models') = list(two_models)

    return (  preds_n )
  } else {
    preds_nmin1 = recursive_ice_long(Tt,
                                     n_nested - 1,
                                     df_obs,
                                     df_interv,
                                     inside_formula_t,
                                     inside_formula_tmin1,
                                     outside_formula,
                                     inside_family,
                                     t_col,
                                     binomial_n,
                                     models,
                                     tmle,
                                     weights,
                                     suppress_rank_warnings)

    obs_keep = t_col >= n_nested #we can only go outward n_nested expectations if that doesn't take us before time 0
    formula_n = glue::glue(outside_formula)

    two_models = list()
    two_models[['tmin1']] = fit_outer_model_long(n_nested, df_obs, preds_nmin1[,1], formula_n, binomial_n, obs_keep)
    two_models[['t']] = fit_outer_model_long(n_nested, df_obs, preds_nmin1[,2], formula_n, binomial_n, obs_keep)

    preds_n = preds_nmin1 #carry outward for observations where t < n_nested
    preds_n[obs_keep, ] = sapply(X = two_models,
                                 FUN = if(suppress_rank_warnings) safe_predict else predict,
                                 newdata=df_interv[obs_keep, ], simplify=TRUE)

    if(tmle) {
      for(s in c('tmin1', 't')) {
        preds_n[obs_keep, s] = tmle_update(preds  = preds_n[obs_keep, s, drop=TRUE],
                                           y      = preds_nmin1[obs_keep, s, drop=TRUE],
                                           w      = df_obs[obs_keep, glue::glue(weights), drop=TRUE],
                                           family = gaussian,
                                           subset = df_obs[obs_keep, glue::glue('A_lag{n}')]==0)
      }
    }

    if(models) {
      attr(preds_n, 'models') = c(attr(preds_n, 'models'), list(two_models))
      names(attr(preds_n, 'models')) = paste0('n_nested', 0:n_nested)
    }

    return ( preds_n )
  }
}




fit_inner_models_long <- function(df_obs, formula_t, formula_tmin1, family) {
  list(tmin1 = glm(formula_tmin1, family, df_obs),
       t     = glm(formula_t,     family, df_obs))

}

fit_outer_model_long <- function(n_nested, df_obs, preds, formula_n, binomial_n=NULL,
                                     obs_keep=NULL #this variable is used to select observations where t >= k, so we don't get into pre-time-zero lags. Is there a way to automate this?
) {

  if(is.null(binomial_n)) {
    df_obs$lm_weights = 1 #these are 1's unless binomial aggregate data
  } else {
    df_obs$lm_weights = binomial_n / sum(binomial_n)
  }

  df_obs$preds = preds

  return (
    lm( formula = formula_n,
        data=df_obs[obs_keep, ],
        weights = lm_weights)
  )
}
estimate_ice_long <- function(ice_preds, #(NxTt)x2 matrix
                         t_col, #(NxTt)x1 vector indicating which times the ice_preds correspond to
                         link_fun = NULL,
                         binomial_n = NULL) {
  if(is.null(link_fun)) {
    link_fun = function(x) x
  } else {
    stopifnot(is.function(link_fun))
  }
  if(is.null(binomial_n)) binomial_n = rep(1, nrow(ice_preds))
  freq_w = binomial_n / sum(binomial_n)

  df = tibble::tibble(t = t_col, ice_t = ice_preds[,2], ice_s = ice_preds[,1], freq_w=freq_w) %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(estimate = link_fun(sum(ice_t * freq_w) / sum(freq_w)) - link_fun(sum(ice_s * freq_w) / sum(freq_w)))

  return ( df$estimate )
}




