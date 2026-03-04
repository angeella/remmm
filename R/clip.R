#' Clusterwise sign-flipping score test
#'
#'
#' Implements a *clusterwise sign-flipping score-based test* for one or multiple outcomes
#' under within-cluster dependence (i.e., non-independent observations), as in
#' linear mixed-model–like settings.
#'
#' @param formula A model \code{\link{formula}}. For multivariate outcomes, \code{formula}
#'   may refer to a matrix response in \code{data} (as in \code{\link[stats]{lm}})
#'   or to a column in \code{data} handled by \code{make_formula_list()}.
#' @param data A \code{data.frame} containing the variables referenced in \code{formula}.
#' @param n_flips Integer. Number of flips/permutations to generate.
#' @param alternative Character string specifying the alternative hypothesis.
#'   Typically one of \code{"two.sided"}, \code{"greater"}, \code{"less"}.
#' @param id Cluster identifier. Either (i) a vector of cluster ids of length \code{nrow(data)},
#'   or (ii) a single character string naming the clustering column in \code{data}.
#' @param seed Optional integer seed for reproducibility. If not \code{NULL},
#'   \code{set.seed(seed)} is called at the beginning of the function.
#' @param covhat Covariance object passed to \code{whiten_mf_from_V()} as \code{V}.
#'   It can be a single matrix (univariate case) or a list of matrices (one per outcome),
#'   depending on the expected input of \code{whiten_mf_from_V()}.
#'
#' @details
#' The function first expands \code{formula} into a list of equations via
#' \code{make_formula_list()}. Flips are generated with \code{flipscores:::.make_flips()}
#' using the supplied cluster identifier \code{id}. For each equation, the model frame
#' is transformed/whitened via \code{whiten_mf_from_V()} using \code{covhat}, then
#' \code{flipscores()} is run on the transformed data after removing the intercept
#' (\code{. ~ . - 1}). Per-equation objects are stored in \code{mods} and aggregated
#' summaries are returned in \code{Tspace} and \code{summary_table}.
#'
#' @return An object of class \code{"jointest"} (a list) with components:
#' \describe{
#'   \item{Tspace}{Aggregated test-statistic space computed from all equations.}
#'   \item{summary_table}{Aggregated summary table across equations.}
#'   \item{mods}{Named list of per-equation \code{flipscores} objects (names are response variables).}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso \code{\link[flipscores]{flipscores}}
#'
#' @import flipscores
#' @import stats
#' @import jointest
#'
#' @examples
#' library(remmm)
#' df <- data.frame(y = rnorm(100), x = rnorm(100), id = rep(1:20, each = 5))
#' V  <- diag(100)
#' out <- clip(y ~ x, data = df, n_flips = 999, alternative = "two.sided",
#'             id = "id", seed = 1, covhat = V)
#' out$summary_table
#'
#'
#' @export
clip <- function(formula, data, n_flips = 5000, alternative= "two.sided", id, seed = NULL, covhat = NULL){

  if(!is.null(seed)) set.seed(seed)

  res <- make_formula_list(formula = formula, data = data)

  df <- res$data
  formula_list <- res$formula_list

  n_obs <- nrow(df)

  m <- length(res$formula_list)
  var_names <- attr(terms(formula_list[[1]]), "term.labels")
  p <- length(var_names)

  if(is.character(id)){
    id <- df[[id]]
  }

  flips <- .make_flips(n_obs = n_obs, n_flips = n_flips, id = id)

  y_names <- vapply(formula_list, function(f) as.character(f[[2]]), character(1))

  mods <- lapply(seq_len(m), function(x){

    unique_ids <- unique(id)

    out <- whiten_mf_from_V(formula = formula_list[[x]],
                            data = data, V = covhat, id_col = "id")

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
