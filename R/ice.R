#' Iterated conditional expectation (ICE) estimator
#'
#' @param df_obs  Data frame with one row per individual-period (if `long=TRUE`) or one row per individual (if `long=FALSE`)
#' @param df_interv Data frame with same dimensions as `df_obs`, but with exposure variables set to intervened values
#' @param inside_formula_t chr, right-hand-side formula for inside model for Yt
#' @param inside_formula_tmin1 chr, right-hand-side formula for inside model for Yt-1
#' @param outside_formula chr, right-hand-side formula for outside models
#' @param Tt int. max periods
#' @param t_col integer vector of length equal to `nrow(df_obs).` Column denoting periods; only needed if `long=TRUE`.
#' @param n_nested int. How many nested expectations should be estimated, starting from the innermost (=0) to the outermost (=Tt)?
#' @param inside_family stats::family object or string referring to one, as in `glm`.
#' @param binomial_n int length nrow(data). Group sizes for binomial aggregate data.
#' @param pt_link_fun function. The scale on which parallel trends is assumed (e.g., `qlogis` for logit scale). Default `NULL` for untransformed scale.
#' @param tibble logical. return results as a tibble (TRUE) or vector (FALSE)?
#' @param models lgl. Return all models as an attribute?
#' @param long lgl. Is df_obs wide (`FALSE`, default) or long (`TRUE`) format?
#'
#'
#' @return
#' @export
#'
#' @examples
ice <- function(df_obs,
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
                long=FALSE,
                suppress_rank_warnings=FALSE) {

  list_args = c(as.list(environment()))
  list_args$long = NULL
  if(!long) list_args$t_col <- list_args$n_nested <- NULL

  if(long) {
    do.call(ice_pipeline_long, list_args)
  } else {
    do.call(ice_pipeline_wide, list_args)
  }
}
