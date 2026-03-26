#' @title Cluster wise sign-flipping score test
#' @description This function implements a *cluster wise sign-flipping score-based test* for
#' one or multiple outcomes under within-cluster dependence (i.e., non-independent observations).
#' @usage clip(formula, data, cluster,
#' n_flips = 5000, alternative= "two.sided", seed = NULL, V = NULL)
#' @param formula A model \code{\link{formula}}. For multivariate outcomes, \code{formula}
#'   may refer to a matrix response in \code{data} (as in \code{\link[stats]{lm}})
#'   or to a column in \code{data}.
#' @param data A \code{data.frame} containing the variables referenced in \code{formula}.
#' @param cluster Cluster identifier. Character string naming the clustering column in \code{data}.
#' @param n_flips Integer. Number of flips to generate. Default is 5000.
#' @param alternative Character string specifying the alternative hypothesis used to test the fixed effect
#' coefficients. One of \code{"two.sided"}, \code{"greater"}, \code{"less"}.
#' Default is \code{"two.sided"}.
#' @param seed Optional integer seed for reproducibility.
#' @param V Working covariance matrix to be used inside the score contributions.
#' It can be a single matrix (univariate case, one response) or a list of matrices (one per outcome).
#' Default is NULL, i.e., identity matrix.
#' @return A \code{jointest} object, i.e., a list containing the following objects:
#' \describe{
#'   \item{Tspace}{\code{data.frame} where rows represents the sign-flipping transformed (plus the identity one) test and columns the variables.}
#'   \item{summary_table}{\code{data.frame} containing for each model the estimated parameter(s), score(s), std error(s), test(s), partial correlation(s) and p-value(s).}
#'   \item{mods}{List of \code{flipscores} objects, one for each outcome.}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso \code{\link[flipscores]{flipscores}}
#'
#' @import flipscores
#' @import stats
#' @import jointest
#' @author Angela Andreella
#' @examples
#' library(remmm)
#' df <- data.frame(y = rnorm(100), x = rnorm(100), id = rep(1:20, each = 5))
#' V  <- diag(100)
#' out <- clip(y ~ x, data = df, cluster = "id", V = V)
#' summary(out)
#'
#'
#' @export
#'
clip <- function(formula, data, cluster, n_flips = 5000, alternative= "two.sided", seed = NULL, V = NULL){

  if(!is.null(seed)) set.seed(seed)

  res <- make_formula_list(formula = formula, data = data)

  df <- res$data
  formula_list <- res$formula_list

  n_obs <- nrow(df)

  m <- length(res$formula_list)
  var_names <- attr(terms(formula_list[[1]]), "term.labels")
  p <- length(var_names)

  cluster_col <- df[[cluster]]


  flips <- .make_flips(n_obs = n_obs, n_flips = n_flips, id = cluster_col)

  y_names <- vapply(formula_list, function(f) as.character(f[[2]]), character(1))

  mods <- lapply(seq_len(m), function(x){

    unique_ids <- unique(cluster_col)

    out <- whiten_mf_from_V(formula = formula_list[[x]],
                            data = data, V = V, id_col = cluster)

    data_star <- as.data.frame(out$mf_star)

    formula_up <- update(formula_list[[x]], . ~ . - 1)

    temp <- flipscores(
      formula = formula_up,
      data    = data_star,
      score_type = "standardized",
      flips = flips,
      alternative = alternative
    )
    temp$summary_table <- .get_summary_table_from_flipscores(temp)
    temp
  })

  names(mods) <- y_names

  out <- list(
    Tspace = .get_all_Tspace(mods),
    summary_table = .get_all_summary_table(mods),
    mods = mods,
    call = match.call()
  )
  class(out) <- unique(c("jointest", class(out)))
  out
}
