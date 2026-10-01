

prepareNullSpace <- function(Zg, tol = 1e-10) {
  Zg <- as.matrix(Zg)
  N <- nrow(Zg)
  q <- ncol(Zg)

  sv <- svd(Zg, nu = min(N, q), nv = min(N, q))
  scale_tol <- tol * max(1, sv$d[1L])
  K <- sum(sv$d > scale_tol)

  U1 <- sv$u[, seq_len(K), drop = FALSE]
  V <- sv$v[, seq_len(K), drop = FALSE]
  d <- sv$d[seq_len(K)]

  Q <- qr.Q(qr(U1), complete = TRUE)
  U2 <- Q[, seq.int(K + 1L, N), drop = FALSE]

  list(
    U1 = U1,
    U2 = U2,
    V = V,
    d = d,
    K = K,
    nullDim = N - K,
    N = N,
    q = q,
    colnames = colnames(Zg)
  )
}

ordered_pair_covariance <- function(r, group) {
  group <- factor(group)
  n <- as.numeric(tapply(r, group, length))
  s <- as.numeric(tapply(r, group, sum))
  ss <- as.numeric(tapply(r^2, group, sum))

  ok <- n >= 2L
  den <- sum(n[ok] * (n[ok] - 1L))
  if (den <= 0) return(NA_real_)

  num <- sum(s[ok]^2 - ss[ok])
  num / den
}

estimate_vc_mom <- function(db, y_col = "y", x_col = "x",
                            subject_col = "subjects", item_col = "items",
                            eps = 1e-8) {
  required_cols <- c(y_col, x_col, subject_col, item_col)
  missing_cols <- setdiff(required_cols, names(db))

  y <- as.numeric(db[[y_col]])
  x <- as.numeric(db[[x_col]])
  X <- matrix(x, ncol = 1L)

  beta_hat <- as.numeric(qr.coef(qr(X), y))
  if (!is.finite(beta_hat)) beta_hat <- 0
  r <- as.numeric(y - X %*% beta_hat)

  sig_total <- mean(r^2)
  covA <- ordered_pair_covariance(r, db[[subject_col]])
  covB <- ordered_pair_covariance(r, db[[item_col]])

  sigA2 <- if (is.finite(covA)) max(covA, eps) else eps
  sigB2 <- if (is.finite(covB)) max(covB, eps) else eps
  sigE2 <- max(sig_total - sigA2 - sigB2, eps)

  list(
    sigA2 = sigA2,
    sigB2 = sigB2,
    sigE2 = sigE2,
    beta_hat = beta_hat,
    covA = covA,
    covB = covB,
    sig_total = sig_total,
    method = "mom"
  )
}

estimate_vc_lmer <- function(db, y_col = "y", x_col = "x",
                             subject_col = "subjects", item_col = "items",
                             eps = 1e-8) {


  required_cols <- c(y_col, x_col, subject_col, item_col)
  missing_cols <- setdiff(required_cols, names(db))
  tmp <- db
  tmp$y <- as.numeric(tmp[[y_col]])
  tmp$x <- as.numeric(tmp[[x_col]])
  tmp$subjects <- factor(tmp[[subject_col]])
  tmp$items <- factor(tmp[[item_col]])

  mod <- lme4::lmer(y ~ x + (1 | subjects) + (1 | items), data = tmp)
  vc <- as.data.frame(lme4::VarCorr(mod))

  get_vc <- function(grp) {
    out <- vc$vcov[vc$grp == grp & (is.na(vc$var1) | vc$var1 == "(Intercept)")]
    if (length(out) == 0L || !is.finite(out[1L])) eps else max(out[1L], eps)
  }

  list(
    sigA2 = get_vc("subjects"),
    sigB2 = get_vc("items"),
    sigE2 = max(vc$vcov[vc$grp == "Residual"][1L], eps),
    method = "lmer"
  )
}

computeDeltaCol <- function(decomp, lambda_vec) {
  K <- decomp$K
  q <- decomp$q

  V <- decomp$V
  D <- diag(decomp$d, nrow = K, ncol = K)
  VD <- V %*% D
  Pl <- diag(lambda_vec, nrow = q, ncol = q)

  A <- tcrossprod(VD) + Pl
  Delta_col <- diag(K) - t(VD) %*% solve(A, VD)
  (Delta_col + t(Delta_col)) / 2
}




random_rotation <- function(n) {
  Z <- matrix(stats::rnorm(n * n), nrow = n, ncol = n)
  qrz <- qr(Z)
  Q <- qr.Q(qrz)
  R <- qr.R(qrz)
  s <- sign(diag(R))
  s[s == 0] <- 1
  sweep(Q, 2L, s, FUN = "*")
}


compute_null <- function(X2, R2, B,transformation = "signflip") {

  X2 <- as.numeric(X2)
  R2 <- as.matrix(R2)
  n <- length(X2)

  L <- ncol(R2)

  if (transformation == "signflip") {
    signs <- matrix(sample(c(-1, 1), B * n, replace = TRUE), nrow = B, ncol = n)
    signs[1L, ] <- 1
    weighted_signs <- sweep(signs, 2L, X2, FUN = "*")
    Mnull <- weighted_signs %*% R2
  }

  if (transformation == "permutation") {
    perms <- matrix(NA_integer_, nrow = B, ncol = n)
    perms[(2:B), ] <- replicate(B-1, sample.int(n))
    perms[1L, ] <- seq_len(n)
    Mnull <- matrix(NA_real_, nrow = B, ncol = L)
    for (b in seq_len(B)) {
      Mnull[b, ] <- as.numeric(crossprod(X2, R2[perms[b, ], , drop = FALSE]))
    }
  }
  if (transformation == "rotation") {
    Mnull <- matrix(NA_real_, nrow = B, ncol = L)
    Mnull[1L, ] <- as.numeric(crossprod(X2, R2))
    if (B >= 2L) {
      for (b in 2L:B) {
        Q <- random_rotation(n)
        Mnull[b, ] <- as.numeric(crossprod(X2, Q %*% R2))
      }
    }
  }
  return(Mnull)

}
