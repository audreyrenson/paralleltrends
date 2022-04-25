#' Targeted maximum likeilhood estimators (TMLE)
#'
#' Functions to estimate the parallel trends g-formula using double robust targeted maximum likelihood estimators.
#'
#' @param df_obs  Data frame with one row per individual-period (if `long=TRUE`) or one row per individual (if `long=FALSE`)
#' @param df_interv Data frame with same dimensions as `df_obs`, but with exposure variables set to intervened values
#' @param den_formula chr. glue-style formula for denominator model(s)
#' @param inside_formula_t chr, glue-style right-hand-side formula for inside model for Yt
#' @param inside_formula_tmin1 chr, glue-style right-hand-side formula for inside model for Yt-1
#' @param outside_formula chr, glue-style right-hand-side formula for outside models
#' @param Tt int. max periods
#' @param t_col integer vector of length equal to `nrow(df_obs).` Column denoting periods; only needed if `long=TRUE`.
#' @param id chr. name of column in `df_obs` corresponding to unit identifier (only need it `long=TRUE`)
#' @param n_nested int. How many nested expectations should be estimated, starting from the innermost (=0) to the outermost (=Tt)?
#' @param den_family stats::family object or string referring to one, as in `glm`, for denominator model(s).
#' @param inside_family stats::family object or string referring to one, as in `glm`, for innnermost model(s).
#' @param binomial_n int length nrow(data). Group sizes for binomial aggregate data.
#' @param pt_link_fun function. The scale on which parallel trends is assumed (e.g., `qlogis` for logit scale). Default `NULL` for untransformed (identity) scale.
#' @param tibble logical. return results as a tibble (TRUE) or vector (FALSE)?
#' @param models lgl. Return all models as an attribute?
#' @param long lgl. Is df_obs wide (`FALSE`, default) or long (`TRUE`) format?
#' @param suppress_rank_warnings lgl. Rank deficient models are often expected in this setting. Option to turn off warning 'prediction from a rank-deficient fit may be misleading'
#'
#'
#' @return tibble with `Tt` rows and 2 columns. Column `estimate` contains estimates of counterfactual trends from t-1 to t. I.e., these are
#' estimates of \eqn{g\{E[Y_t(\bar a^*)]\} - g\{E[Y_{t-1}(\bar a^*)]\}}, where \eqn{g\{\cdot\}} is the parallel trends link function specified
#' by `pt_link_fun`.
#' @export
#'
#' @examples
#' Tt = 3
#' N = 100
#' Beta = generate_parameters(Tt=Tt)
#' df_obs = generate_data(N, Tt, Beta)
#' df_interv = df_obs %>% dplyr::mutate(A1=0, A2=0, A3=0)
#' tmle(df_obs=df_obs,
#'     df_interv = df_interv,
#'     den_formula = 'A{t}~A{t-1}*(W1{t}+W2{t}+I(W2{t}^2))',
#'     inside_formula_t = '~A{t}*(W1{t}+W2{t}+I(W2{t}^2))',
#'     inside_formula_tmin1 = '~A{t}*(W1{t-1}+W2{t-1}+I(W2{t-1}^2))',
#'     outside_formula = '~A{k}*(W1{k}+W2{k}+I(W2{k}^2))',
#'     Tt=Tt)
tmle = function(df_obs,
                df_interv,
                den_formula,
                inside_formula_t,
                inside_formula_tmin1,
                outside_formula,
                Tt,
                t_col,
                id,
                n_nested=Tt,
                den_family='binomial',
                inside_family='gaussian',
                pt_link_fun=NULL,
                binomial_n=NULL,
                tibble=TRUE,
                models=TRUE,
                long=FALSE,
                suppress_rank_warnings=FALSE)

{
  list_args = c(as.list(environment()))
  list_args$long = NULL
  if(!long) list_args$t_col <- list_args$id <- list_args$n_nested <- NULL #how to incorporate n_nested into wide format?

  if(long) {
    do.call(tmle_pipeline_long, list_args)
  } else {
    do.call(tmle_pipeline_wide, list_args)
  }

}

tmle_pipeline_wide <- function(df_obs,
                               df_interv,
                               Tt,
                               den_formula,
                               inside_formula_t,
                               inside_formula_tmin1,
                               outside_formula,
                               den_family,
                               inside_family,
                               pt_link_fun=NULL,
                               binomial_n=NULL,
                               tibble=TRUE,
                               models=TRUE,
                               suppress_rank_warnings=FALSE)
{

  list_args = c(as.list(environment()))
  list_args$den_formula <- list_args$den_family <- NULL
  list_args$tmle = TRUE

  #calculate TMLE weights
  den_glued_formulas = lapply(1:Tt, function(k) glue_formula(den_formula,t=k))
  den_models = lapply(den_glued_formulas, stats::glm, family=den_family, data=df_obs)
  den_probs = cbind(1, sapply(den_models,
                              if(suppress_rank_warnings) safe_predict else stats::predict,
                              newdata=df_interv,
                              type='response'))
  den = matrixStats::rowCumprods(den_probs)
  a_vars = sapply(den_glued_formulas, function(f) as.character(f)[2])
  num = cbind(1, matrixStats::rowCumprods(df_obs[, a_vars] == df_interv[, a_vars]))
  list_args$weights = num/den

  #run ICE pipeline with tmle=TRUE
  do.call(ice_pipeline_wide, list_args)



}


tmle_pipeline_long <- function(df_obs,
                               df_interv,
                               den_formula,
                               inside_formula_t,
                               inside_formula_tmin1,
                               outside_formula,
                               Tt,
                               t_col,
                               id,
                               n_nested=Tt,
                               den_family,
                               inside_family,
                               pt_link_fun=NULL,
                               binomial_n=NULL,
                               tibble=TRUE,
                               models=TRUE,
                               suppress_rank_warnings=FALSE) {

  list_args = c(as.list(environment()))
  list_args$den_formula <- list_args$den_family <- list_args$id <- NULL
  list_args$tmle = TRUE

  #calculate TMLE weights
  den_formula_glued = glue_formula(den_formula, tvars=paste(attr(df_obs, 'timevars')[-1], collapse='+'))
  den_mod = stats::glm(den_formula_glued, family=den_family, data=df_obs)
  exposure_variable = stats::terms(den_formula_glued)[[2]]
  den_probs = dens(df_interv[[exposure_variable]], model=den_mod, newdata=df_interv, suppress_rank_warnings)

  list_args$weights = df_obs %>%
    dplyr::mutate(den_t=ifelse(t_col==0, 1, den_probs),
                  num_t=.data[[exposure_variable]] == df_interv[[exposure_variable]]) %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::mutate(weight = cumprod(.data[['num_t']]/.data[['den_t']])) %>%
    `[[`('weight')

  #run ICE pipeline with tmle=TRUE
  do.call(ice_pipeline_long, list_args)

}
tmle_update_wide <- function(preds, y, w, family) {
  df = data.frame(y,preds,w)
  pr = safe_predict(stats::glm(y ~ offset(preds), family=family, weights=w), type='response', newdata=df)
  return ( pr )
}
tmle_update_long <- function(preds, y, df_tvars, w, family) {
  df = data.frame(y,preds,w,df_tvars)
  mod = stats::glm(y ~ .-preds-w + offset(preds), family=family, weights=w, data=df)
  pr = safe_predict(mod, type='response', newdata=df)
  return(pr)
}
