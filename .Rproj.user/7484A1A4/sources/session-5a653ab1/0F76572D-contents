make_formula_list <- function(formula, data) {

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  # nome "pulito" della risposta (serve per Y.1, Y.2, ...)
  lhs_name <- deparse(lhs)

  # prova a valutare LHS dentro data; se fallisce, prova nell'ambiente
  y_val <- tryCatch(
    eval(lhs, envir = data, enclos = parent.frame()),
    error = function(e) eval(lhs, envir = parent.frame())
  )

  if (!is.matrix(y_val)) {
    stop("Il lato sinistro non è una matrice (o non è stato trovato).")
  }

  m <- ncol(y_val)

  # creo Y.1, Y.2, ...
  y_names <- paste0(lhs_name, ".", seq_len(m))

  # aggiungo le colonne al data
  for (j in seq_len(m)) {
    data[[y_names[j]]] <- y_val[, j]
  }

  # RHS come stringa
  rhs_str <- paste(deparse(rhs), collapse = "")

  # formule
  formula_list <- lapply(seq_len(m), function(j) {
    as.formula(paste0(y_names[j], " ~ ", rhs_str))
  })

  list(formula_list = formula_list, data = data)
}
