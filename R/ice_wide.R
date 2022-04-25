ice_pipeline_wide <- function(df_obs,
                              df_interv,
                              inside_formula_t,
                              inside_formula_tmin1,
                              outside_formula,
                              Tt,
                              inside_family='gaussian',
                              pt_link_fun=NULL,
                              binomial_n=NULL,
                              tibble=TRUE,
                              models=TRUE,
                              tmle=FALSE,
                              weights=NULL,
                              suppress_rank_warnings=FALSE) {
  # @param df_obs Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
  # @param df_interv `df_obs` with exposures set to their intervened values.
  # @param inside_formula_t chr, right-hand-side formula for inside model for Yt
  # @param inside_formula_tmin1 chr, right-hand-side formula for inside model for Yt-1
  # @param outside_formula chr, right-hand-side formula for outside models
  # @param Tt int. max periods
  # @param inside_family stats::family object or string referring to one, as in `glm`.
  # @param binomial_n int length nrow(df_obs). Group sizes for binomial aggregate data.
  # @param pt_link_fun function. The scale on which parallel trends is assumed (e.g., `qlogis` for logit scale). Default `NULL` for untransformed scale.
  # @param tibble logical. return results as a tibble (TRUE) or vector (FALSE)?
  # @param models logical. Return all models as an attribute?
  # @param tmle lgl. Incorporate tmle updating step to make estimator doubly robust?
  # @param weights Matrix of IPTW weights as returned by `cbind(1, calc_weights(...))`. Default `NULL`, required if `tmle=TRUE`.
  # @param suppress_rank_warnings Rank deficient models are often expected in this setting. Option to turn off warning 'prediction from a rank-deficient fit may be misleading'

  ice_preds = lapply(1:Tt, function(t) recursive_ice_wide(t, k=0,
                                                          df_obs=df_obs,
                                                          df_interv=df_interv,
                                                          inside_formula_t=inside_formula_t,
                                                          inside_formula_tmin1=inside_formula_tmin1,
                                                          outside_formula=outside_formula,
                                                          inside_family=inside_family,
                                                          binomial_n=binomial_n,
                                                          models=models,
                                                          tmle=tmle,
                                                          weights=weights,
                                                          suppress_rank_warnings=suppress_rank_warnings))

  ice_estimates = sapply(ice_preds, estimate_ice_wide, link_fun=pt_link_fun, binomial_n)

  if(!tibble) {
    result = ice_estimates
  } else {
    result = tibble::tibble(t=1:Tt,
                            estimate = ice_estimates)
  }

  if(models) {
    attr(result, 'models') = lapply(ice_preds, function(x) attr(x, 'models'))
    names(attr(result, 'models')) = paste0('t', 1:Tt)
  }

  return(result)

}

recursive_ice_wide <- function(t, k, df_obs, df_interv, inside_formula_t, inside_formula_tmin1,
                               outside_formula, inside_family, binomial_n=NULL, models=TRUE,
                               tmle=FALSE, weights=NULL, suppress_rank_warnings=TRUE) {
  # Recursively apply the ICE algorithm for wide data
  #
  # @param t int. Timepoint the estimand applies to
  # @param k int. Always 0 if called at the user level. Timepoint the outer expectation applies to
  # @param df_obs wide data frame
  # @param df_interv `df_obs` with exposures set to their intervened values.
  # @param inside_formula_t chr. glue-style formula for inside model for Yt (times indexed by t, e.g. '~L{t}')
  # @param inside_formula_tmin1 chr. glue-style formula for inside model for Yt (times indexed by t, e.g. '~L{t-1}')
  # @param outside_formula chr. glue-style formula for outside models (times indexed by k, e.g. '~L{k})
  # @param inside_family stats::families object (or chr) for glm inside models
  # @param binomial_n vector of length nrow(df_obs) of group sizes if Y is binomial aggregate data
  # @param models lgl. Return all models as an attribute?
  # @param tmle lgl. Incorporate tmle updating step to make estimator doubly robust?
  # @param weights Matrix of IPTW weights as returned by `cbind(1, calc_weights(...))`. Default `NULL`, required if `tmle=TRUE`.


  if(k==t) {
    two_models = fit_inner_models_wide(t=k,
                                       df_obs=df_obs,
                                       rhs_formula_t = inside_formula_t,
                                       rhs_formula_tmin1 = inside_formula_tmin1,
                                       family=inside_family,
                                       binomial_n)


    preds_k = sapply(X=two_models,
                     FUN=if(suppress_rank_warnings) safe_predict else stats::predict,
                     newdata=df_interv,
                     type='response',
                     simplify=TRUE)

    if(tmle) {### BOOKMARK 4.21.22. Should be last step to making tmle_pipeline_wide() work - make variable labels more general and use df_interv to find subset.
      if(is.null(weights)) stop('must supply weights if tmle=TRUE')
      preds_k[,1] = tmle_update_wide(preds_k[,1],
                                     y=df_obs[[glue::glue('Y{t-1}')]],
                                     w=weights[,k+1],
                                     family=inside_family)
      preds_k[,2] = tmle_update_wide(preds_k[,2],
                                     y=df_obs[[glue::glue('Y{t}')]],
                                     w=weights[,k+1],
                                     family=inside_family)
    }

    if(models)  attr(preds_k, 'models') = list(two_models)

    return (  preds_k )
  } else {
    preds_kplus1 = recursive_ice_wide(t,
                                      k+1,
                                      df_obs,
                                      df_interv,
                                      inside_formula_t,
                                      inside_formula_tmin1,
                                      outside_formula,
                                      inside_family,
                                      binomial_n,
                                      models,
                                      tmle,
                                      weights,
                                      suppress_rank_warnings)

    two_models = lapply(1:2, function(h)  fit_outer_model_wide(k, df_obs, preds_kplus1[,h], outside_formula, binomial_n))

    preds_k = sapply(X=two_models,
                     FUN=if(suppress_rank_warnings) safe_predict else stats::predict,
                     newdata=df_interv,
                     type='response',
                     simplify=TRUE)

    if(tmle) {
      preds_k[,1] = tmle_update_wide(preds_k[,1],
                                     y=preds_kplus1[,1],
                                     w=weights[,k+1],
                                     family=inside_family)
      preds_k[,2] = tmle_update_wide(preds_k[,2],
                                     y=preds_kplus1[,2],
                                     w=weights[,k+1],
                                     family=inside_family)
    }


    if(models) {
      names(two_models) = c('tmin1','t')
      attr(preds_k, 'models') = c(attr(preds_kplus1,'models'), list(two_models))
      names(attr(preds_k, 'models')) = paste0('n_nested', 0:(t-k))
    }


    return ( preds_k )
  }
}


fit_inner_models_wide <- function(t, df_obs, rhs_formula_t, rhs_formula_tmin1, family, binomial_n=NULL) {
  #this is for estimating the innermost expectation of ICE,
  #two at once because we typically want E[Y_t|L_t, A_t] and E[Y_{t-1}|L_t, A_t]

  # @param t int.
  # @param df_obs wide data frame
  # @param rhs_formula_t chr. RHS glue-style formula for Yt model
  # @param rhs_formula_tmin chr. RHS glue-style formula for Yt-1 model
  # @param family stats::families object or string referring to one
  # @param binomial_n vector of group sizes of length=nrow(df_obs)


  if(is.null(binomial_n)) { #y is not aggregate binomial
    formulas =     c(glue_formula(paste0('Y{t-1}', rhs_formula_tmin1), t=t),
                     glue_formula(paste('Y{t}',rhs_formula_t ), t=t))
  } else {#y is aggregate binomial
    formulas =     c(glue_formula(paste0('cbind(Y{t-1}, binomial_n - Y{t-1})', rhs_formula_tmin1), t=t),
                     glue_formula(paste('cbind(Y{t}, binomial_n - Y{t})',rhs_formula_t ), t=t))
    df_obs$binomial_n = binomial_n #and we don't want any ambiguity in case the 'df_obs' already has a column 'binomial_n' (!!)
  }

  models = lapply(1:2, function(s) stats::glm(formula = stats::as.formula(formulas[[s]]),
                                              family,
                                              data=df_obs))

  names(models) = c('tmin1', 't')

  return (models)
}

fit_outer_model_wide <- function(k, df_obs, preds, rhs_formula, binomial_n=NULL) {
  # Fit outer models for ICE using `lm`
  #
  # @param k int. Timepoint to which the outer expectation, i.e. E(mu_k^t|L_k, A_k), applies.
  # @param df_obs wide data frame
  # @param preds vector of length(nrow(df_obs)) of values predicted by the k+1 model.
  # @param rhs_formula_t chr. RHS glue-style formula for Yt model
  # @param rhs_formula_tmin chr. RHS glue-style formula for Yt-1 model
  # @param binomial_n vector of group sizes if Y refers to binomial aggregate data. Used to construct weights.


  #t is the timepoint for which we are estimating \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)]
  #k is the the timepoint the outer expectation applies to
  #preds are the k+1th outer expectation

  if(is.null(binomial_n)) {
    lm_weights = 1 #these are 1's unless binomial aggregate data
  } else {
    lm_weights = binomial_n / sum(binomial_n)
  }

  df_obs$preds = preds
  df_obs$lm_weights = lm_weights

  mod = stats::lm( formula = glue_formula(paste0('preds', rhs_formula), k=k),
                   data=df_obs,
                   weights = lm_weights)

  return (mod)
}


estimate_ice_wide <- function(ice_preds,  # Nx2 matrix
                         link_fun,
                         binomial_n = NULL) {

  if(is.null(link_fun)) {
    link_fun = function(x) x
  } else {
    stopifnot(is.function(link_fun))
  }

  if(is.null(binomial_n)) binomial_n = rep(1, nrow(ice_preds))
  freq_w = binomial_n / sum(binomial_n)

  return ( diff(link_fun(colSums(ice_preds * freq_w) )) )
}


