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

  model.frame(formula, data = data, ...)
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
