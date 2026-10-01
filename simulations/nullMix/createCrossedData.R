createCrossedData <- function(sA = 2, sB = 2, sE = 1,
                              R, C,
                              beta,
                              seed = 1234,
                              rho = 0) {

  set.seed(seed)

  item_pool <- max(C)
  rows <- vector("list", R)

  for (i in seq_len(R)) {
    items_i <- sample(seq_len(item_pool), size = C[i], replace = FALSE)
    rows[[i]] <- data.frame(subjects = i,items = items_i)
  }

  db <- do.call(rbind, rows)

  db$subjects <- factor(db$subjects, levels = seq_len(R))
  db$items <- factor(db$items, levels = seq_len(item_pool))

  rownames(db) <- NULL

  ZA <- stats::model.matrix( ~ 0 + subjects, data = db)
  ZB <- stats::model.matrix( ~ 0 + items, data = db)
  Zg <- cbind(ZA, ZB)

  N <- nrow(db)
  L <- length(beta)

  covX <- diag(N)
  covX <- (covX + t(covX)) / 2

  x <- t(chol(covX)) %*% rnorm(N)
  db$x <- as.numeric((x - mean(x)) / sd(x))

  if (L == 1) {

    SigmaA <- matrix(sA^2, nrow = 1L, ncol = 1L)
    SigmaB <- matrix(sB^2, nrow = 1L, ncol = 1L)
    SigmaE <- matrix(sE^2, nrow = 1L, ncol = 1L)

  } else {

    lower_bound <- -1 / (L - 1)
    Corr <- matrix(rho,nrow = L,ncol = L)
    diag(Corr) <- 1

    SigmaA <- diag(sA, nrow = L, ncol = L) %*% Corr %*% diag(sA, nrow = L, ncol = L)

    SigmaB <- diag(sB, nrow = L, ncol = L) %*% Corr %*% diag(sB, nrow = L, ncol = L)

    SigmaE <- diag(sE, nrow = L, ncol = L) %*% Corr %*% diag(sE, nrow = L, ncol = L)
  }

  a_effect <- matrix(rnorm(R * L), nrow = R, ncol = L) %*% chol(SigmaA)
  b_effect <- matrix(rnorm(item_pool * L), nrow = item_pool, ncol = L) %*% chol(SigmaB)

  subj_id <- as.integer(db$subjects)
  item_id <- as.integer(db$items)

  eps <- matrix(rnorm(N * L), nrow = N, ncol = L) %*% chol(SigmaE)

  mu <- a_effect[subj_id, , drop = FALSE] + b_effect[item_id, , drop = FALSE]

  signal <- sweep(matrix(db$x, nrow = N, ncol = L), MARGIN = 2, STATS = beta, FUN = "*")

  Y <- mu + signal + eps

  y_cols <- paste0("y", seq_len(L))
  colnames(Y) <- y_cols

  for (ell in seq_len(L)) {
    db[[y_cols[ell]]] <- as.numeric(Y[, ell])
  }

  if (L == 1) {db$y <- as.numeric(Y[, 1])}

  list(db = db, Zg = Zg)
}
