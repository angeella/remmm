
source("simulations/nullMix/utils.R")

nullMix <- function(db,
                    B = 1000,
                    Zg = NULL,
                    y_cols = NULL,
                    returnNullDistr = FALSE,
                    beta0 = 0,
                    known = FALSE,
                    sA = NULL, sB = NULL, sE = NULL,
                    variance_method = "mom",
                    transformation = "signflip") {



  Y <- as.matrix(db[, y_cols])
  X <- db[, "x"]
  N <- length(X)
  L <- length(y_cols)

  decomp <- prepareNullSpace(Zg)

  X1 <- as.numeric(crossprod(decomp$U1, X))
  X2 <- as.numeric(crossprod(decomp$U2, X))

  t_col <- numeric(L)
  R2 <- matrix(NA_real_, nrow = decomp$nullDim, ncol = L, dimnames = list(NULL, y_cols))



  for (ell in seq_len(L)) {
    y <- as.numeric(Y[, ell])

    if (!known) {
      tmp <- db
      tmp$.eigen_y <- y
      if (variance_method == "mom") {
        vc <- estimate_vc_mom(
          tmp,
          y_col = ".eigen_y")
      } else {
        vc <- estimate_vc_lmer(
          tmp,
          y_col = ".eigen_y")
      }
      sigA2 <- vc$sigA2
      sigB2 <- vc$sigB2
      sigE2 <- vc$sigE2
    }else{
      sigA2 <- sA^2
      sigB2 <- sB^2
      sigE2 <- sE^2
    }

    lambdaA <- sigE2 / sigA2
    lambdaB <- sigE2 / sigB2

    cn <- colnames(Zg)
    is_subject <- startsWith(cn, "subjects")

    lambda_vec <- ifelse(is_subject, lambdaA, lambdaB)
    Delta_col <- computeDeltaCol(decomp, lambda_vec)

    r0 <- as.numeric(y - beta0[ell] * X)
    r1 <- as.numeric(crossprod(decomp$U1, r0))
    r2 <- as.numeric(crossprod(decomp$U2, r0))

    t_col[ell] <- as.numeric(crossprod(X1, Delta_col %*% r1))
    R2[, ell] <- r2
  }

  null_part <- compute_null(
    X2 = X2,
    R2 = R2,
    B = B,
    transformation = transformation
  )

  Tmat <- sweep(null_part, 2L, t_col, FUN = "+")

  if(L == 1){
    pv <- mean(abs(Tmat) >= abs(Tmat[1]))
  }else{
    pv <- jointest:::maxT.light(abs(Tmat))
  }
  return(list(pv = pv))

}
