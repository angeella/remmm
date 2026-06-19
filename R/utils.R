make_formula_list <- function(formula, data) {

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  lhs_name <- deparse(lhs)

  y_val <- tryCatch(
    eval(lhs, envir = data, enclos = parent.frame()),
    error = function(e) eval(lhs, envir = parent.frame())
  )


  m <- ncol(y_val)
  if(is.null(m)){

    formula_list <- list(formula)
    data <- data

  }else{
  y_names <- paste0(lhs_name, ".", seq_len(m))

  for (j in seq_len(m)) {
    data[[y_names[j]]] <- y_val[, j]
  }

  rhs_str <- paste(deparse(rhs), collapse = "")

  formula_list <- lapply(seq_len(m), function(j) {
    as.formula(paste0(y_names[j], " ~ ", rhs_str))
  })
  }
  list(formula_list = formula_list, data = data)
}


safe_model_frame <- function(formula, data, ...) {
  vars <- all.vars(formula)

  for (v in setdiff(vars, names(data))) {
    m <- regexec("^(.+)\\.(\\d+)$", v)
    parts <- regmatches(v, m)[[1]]
    if (length(parts) == 3) {
      base <- parts[2]
      k <- as.integer(parts[3])

      if (base %in% names(data)) {
        obj <- data[[base]]

        if (is.matrix(obj) || is.data.frame(obj)) {
          if (ncol(obj) >= k) {
            data[[v]] <- obj[, k]
          }
        }
      }
    }
  }

  data_safe <-cbind(data[, vars[1]], model.matrix(formula, data = data))
  colnames(data_safe)[1] <- vars[1]
  data_safe
}


inv_sqrt_spd <- function(V, jitter = 1e-10, tol_mult = 1e-10) {
  V <- as.matrix(V)
  n <- nrow(V)
  if (n != ncol(V)) stop("V must be square")

  V <- V + jitter * diag(n)

  eig <- eigen(V, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  tol <- tol_mult * max(1, max(vals))
  vals[vals < tol] <- tol

  vecs %*% diag(1 / sqrt(vals), n) %*% t(vecs)
}

left_multiply_blocks <- function(X, id, blocks) {
  X <- if (is.vector(X)) matrix(X, ncol = 1) else as.matrix(X)
  res <- matrix(0, nrow = nrow(X), ncol = ncol(X))

  u <- unique(id)
  for (subj in u) {
    idx <- which(id == subj)
    Vi  <- blocks[[as.character(subj)]]
    res[idx, ] <- Vi %*% X[idx, , drop = FALSE]
  }

  if (ncol(res) == 1) as.numeric(res) else {
    colnames(res) <- colnames(X)
    res
  }
}

is_block_diagonal_by_id <- function(V, id, rel_tol = 1e-12) {
  V <- as.matrix(V)
  n <- nrow(V)
  if (n != ncol(V)) stop("V must be square")
  if (length(id) != n) stop("id length must match nrow(V)")

  # scala “tipica”
  scale <- max(1, max(abs(diag(V))))
  mask_off <- outer(id, id, FUN = "!=")

  max_off <- if (any(mask_off)) max(abs(V[mask_off])) else 0
  max_off <= rel_tol * scale
}

whiten_mf_from_V <- function(formula, data, V, id_col = NULL,
                             rel_tol_block = 1e-12,
                             jitter = 1e-10, tol_mult = 1e-10) {
  mf <- safe_model_frame(formula, data)
  n  <- nrow(mf)

  V <- as.matrix(V)


  use_blocks <- FALSE
  blocks <- NULL
  Vinvsqrt <- NULL

  if (!is.null(id_col)) {
    id <- data[[id_col]]

    use_blocks <- is_block_diagonal_by_id(V, id, rel_tol = rel_tol_block)

    if (use_blocks) {
      u <- unique(id)
      blocks <- setNames(vector("list", length(u)), as.character(u))
      for (subj in u) {
        idx <- which(id == subj)
        V_sub <- V[idx, idx, drop = FALSE]
        blocks[[as.character(subj)]] <- inv_sqrt_spd(V_sub, jitter = jitter, tol_mult = tol_mult)
      }
      mf_star <- left_multiply_blocks(mf, id = id, blocks = blocks)
      return(list(
        mf = mf,
        mf_star = mf_star,
        mode = "block",
        blocks = blocks
      ))
    }
  }

  Vinvsqrt <- inv_sqrt_spd(V, jitter = jitter, tol_mult = tol_mult)
  mf_star <- as.matrix(Vinvsqrt %*% as.matrix(mf))
  colnames(mf_star) <- colnames(mf)

  list(
    mf = mf,
    mf_star = mf_star,
    mode = "global",
    Vinvsqrt = Vinvsqrt
  )
}



.make_summary_table <- function(scores,Tspace,alternative){
  p.values=apply(Tspace,2,flipscores:::.t2p_only_first, alternative)
  data.frame(response=colnames(scores),score= colSums(scores),p.values=p.values)
}


.make_output_from_list_Tspace_summary_table <- function(res_list,original_call){
  Tspace=lapply(res_list,function(x) x$Tspace)
  summary_table=lapply(res_list,function(x) x$summary_table)
  mods=lapply(res_list,function(x) x$mod)
  names(mods)=paste0("mod",1:length(mods))

  Tspace=do.call(cbind,Tspace)
  summary_table=do.call(rbind,summary_table)
  out=list(Tspace=Tspace,
           summary_table=summary_table,
           mods=mods,
           call = original_call)

  class(out) <- c("remmm", class(out))
  class(out) <- c("joint_flipscores", class(out))
  return(out)
}

################

.get_IH <- function (Z)
{
  diag(nrow(Z)) - .get_H(Z)
}
.get_H <- function (Z)
{
  Z %*% solve(t(Z) %*% Z) %*% t(Z)
}

#########
#######################################
#' @examples
#' # example code
#'
#' # scores has only 3 rows
#' scores <- matrix(
#'   c(1, 2, 3, 4, 5, 6, 7, 8, 9),
#'     nrow = 3,
#'       dimnames = list(c("gene1", "gene3", "gene5"), c("S1", "S2", "S3"))
#'       )
#'       scores
#'       #       S1 S2 S3
#'       # gene1  1  4  7
#'       # gene3  2  5  8
#'       # gene5  3  6  9
#'       # cluster_names has 5 elements — more than nrow(scores)
#'       cluster_names <- c("gene1", "gene2", "gene3", "gene4", "gene5")
#'
#' result <- fill_scores_by_cluster(scores, cluster_names)
#' result
#' #       S1 S2 S3
#' # gene1  1  4  7   # copied from scores
#' # gene2  0  0  0   # not in scores → zeroed
#' # gene3  2  5  8   # copied from scores
#' # gene4  0  0  0   # not in scores → zeroed
#' # gene5  3  6  9   # copied from scores
#'
#' cluster_names <- c("gene1", "gene3", "gene5")
#'
#' result <- fill_scores_by_cluster(scores, cluster_names)
#' result
#'
#'@noRd
#'@keywords internal
fill_scores_by_cluster <- function(scores_A, cluster_names) {

  # --- Input Validation ---
  if (!is.matrix(scores_A$scores))
    stop("'scores' must be a matrix.")
  if (is.null(rownames(scores_A$scores)))
    stop("'scores' matrix must have rownames.")
  if (!is.character(cluster_names))
    cluster_names=as.character(cluster_names)

  # Warn if some cluster_names are not found in rownames(scores)
  missing_names <- cluster_names[!cluster_names %in% rownames(scores_A$scores)]
  if (length(missing_names) > 0){
    # --- Build output matrix filled with zeros ---
    # Rows = all cluster_names, Cols = same as scores
    result_scores <- matrix(
      0,
      nrow     = length(cluster_names),
      ncol     = ncol(scores_A$scores),
      dimnames = list(cluster_names, colnames(scores_A$scores))
    )
    # --- Fill in rows that exist in scores ---
    # Only copy rows whose names appear in cluster_names
    matching_names <- cluster_names[cluster_names %in% rownames(scores_A$scores)]
    result_scores[matching_names, ] <- scores_A$scores[matching_names, , drop = FALSE]

    result_A <- matrix(
      0,
      nrow     = length(cluster_names),
      ncol     = ncol(scores_A$A),
      dimnames = list(cluster_names, colnames(scores_A$A))
    )
    result_A[matching_names, ] <- scores_A$A[matching_names, , drop = FALSE]

    return(list(result_scores,result_A))
  } else {
    scores_A$scores=scores_A$scores[cluster_names,,drop=FALSE]
    scores_A$A=scores_A$A[cluster_names,,drop=FALSE]
    return(scores_A)
  }
}



##################
#' @examples
#' # --- Example usage ---
#' D <- data.frame(y1 = c(1, 2, 3, 4),
#'                 y2 = 4:1,
#'                 x = c(2, 4, 6, 8),
#'                 z = c(1, 0, 1, 0))
#'
#' result <- formula_to_matrices(y1 ~ x + z, data = D)
#' result
#'
#' result <- formula_to_matrices(cbind(y1,y2) ~ x * z, data = D)
#' result
#'
#' result$Y  # response matrix
#' result$X  # design matrix (x, z + intercept column)
#'@noRd
#'@keywords internal

formula_to_matrices <- function(formula, data) {

  # Build the model frame (handles NA, subset, etc.)
  mf <- model.frame(formula, data = data)

  # Right-hand side: design matrix (X), includes intercept by default
  X <- model.matrix(formula, data = mf)

  # Left-hand side: response matrix (Y)
  Y <- model.response(mf)
  if(is.vector(Y)) {
    Y= as.matrix(Y)
    colnames(Y)=as.character(formula[[2]])
  }

  list(Y = Y, X = X)
}



