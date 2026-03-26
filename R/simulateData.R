#' @title Simulate multivariate hierarchical data
#'
#' @description
#' Generates multivariate hierarchical data with \eqn{m} outcomes, cluster-level
#' random effects, and correlated residual errors across outcomes.
#'
#' For each outcome \eqn{i = 1, \dots, m}, the response is generated as:
#'
#' \deqn{
#' y_i = X \beta_i + Z \gamma + U + RX \cdot X + RZ \cdot Z + \varepsilon_i
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{X} and \eqn{Z} are fixed-effect covariates,
#'   \item \eqn{U}, \eqn{RX}, \eqn{RZ} are cluster-level random effects,
#'   \item \eqn{\varepsilon} follows a multivariate normal distribution
#'         with covariance matrix \code{Sigma}.
#' }
#'
#' @param Sigma Either:
#'   \itemize{
#'     \item \code{"autocorrelation"}: AR(1)-like structure with correlation \code{rho},
#'     \item \code{"equicorrelation"}: constant correlation \code{rho},
#'     \item \code{"random"}: random positive-definite matrix,
#'     \item or a numeric covariance matrix.
#'   }
#' @param rho Correlation parameter used when \code{Sigma} is
#'   \code{"autocorrelation"} or \code{"equicorrelation"}.
#' @param beta Numeric vector of length \eqn{m}; outcome-specific slopes for \code{X}.
#' @param gamma Numeric scalar; common slope for \code{Z} across outcomes.
#' @param J Integer; number of clusters.
#' @param nJ Integer vector of length \code{J}; number of observations per cluster.
#' @param sd_within Numeric; standard deviation of cluster-level random effects
#' @param seed Integer; random seed for reproducibility.
#' @param sd_eps Numeric; marginal standard deviation used to scale the diagonal
#'   of the residual covariance matrix \code{Sigma}.
#' @param shared_re Logical; if \code{TRUE}, the same random effects are shared
#'   across outcomes; otherwise outcome-specific random effects are generated.
#' @param shared_fe Logical; if \code{TRUE}, covariates \code{X} and \code{Z}
#'   are shared across outcomes; otherwise outcome-specific covariates are generated.
#'
#' @return A list of length \eqn{m}. Each element is a \code{data.frame} with columns:
#' \describe{
#'   \item{\code{id}}{Cluster identifier}
#'   \item{\code{X}, \code{Z}}{Fixed-effect covariates}
#'   \item{\code{U}, \code{RX}, \code{RZ}}{Random effects}
#'   \item{\code{eps}}{Residual term}
#'   \item{\code{y}}{Generated outcome}
#' }
#'
#' The returned list contains two attributes:
#' \describe{
#'   \item{\code{Sigma_eps_true}}{True residual covariance matrix}
#'   \item{\code{Sigma_random_true}}{True random-effects covariance matrix}
#' }
#'
#' @examples
#' \dontrun{
#' library(MASS)
#'
#' db <- simulateData(
#'   Sigma = "equicorrelation",
#'   rho = 0.4,
#'   beta = c(0.5, 1.2),
#'   gamma = 0.8,
#'   J = 8,
#'   nJ = rep(25, 8),
#'   sd_within = 1,
#'   seed = 123,
#'   sd_eps = 1,
#'   shared_re = TRUE,
#'   shared_fe = TRUE
#' )
#'
#' str(db[[1]])
#' attributes(db)
#' }
#'
#' @importFrom MASS mvrnorm
#' @export


simulateData <- function(Sigma, rho,
                         beta, gamma,
                         J, nJ, sd_within,
                         seed, sd_eps,
                         shared_re = TRUE, shared_fe = TRUE){

  set.seed(seed)

  m <- length(beta)

  if(Sigma == "autocorrelation"){
    Sigma <- outer(1:m, 1:m, function(i, j) rho^abs(i - j))
    diag(Sigma) <- sd_eps
  }else{
    if(Sigma == "equicorrelation"){
      Sigma <- matrix(rho, m, m)
      diag(Sigma) <- sd_eps^2
    }else{

      if(Sigma == "random"){
        A <- matrix(runif(m^2)*2-1, ncol=m)
        Sigma <- t(A) %*% A
      }else{
        Sigma <- diag(m)
      }
    }
  }

  db <- array(NA, dim = c(sum(nJ), 7, m))

  db[,1,] <- rep(seq(J), times = nJ)


  if(shared_fe){
    Sigma_fixed <- matrix(c(1, 0.7, 0.7, 1), ncol = 2)

    fixed_effects <- mvrnorm(nrow(db[,,1]), rep(0, 2),
                             Sigma = Sigma_fixed)

    db[,2,] <- matrix(rep(fixed_effects[,1],m), ncol = m)
    db[,3,] <- matrix(rep(fixed_effects[,2],m), ncol = m)
    db[,7,] <- mvrnorm(nrow(db), c(rep(0, m)), Sigma)

  }else{
    Sigma_fixed <- matrix(0.7, 2*m, 2*m)
    diag(Sigma_fixed) <- 1

    fixed_effects <- mvrnorm(nrow(db[,,1]), rep(0, 2*m),
                             Sigma = Sigma_fixed)

    db[,2,] <- fixed_effects[,1:m]
    db[,3,] <- fixed_effects[,(m+1):(2*m)]

    db[,7,] <- mvrnorm(nrow(db), c(rep(0, m)), Sigma)

  }

  Sigma_random <- matrix(c(sd_within^2, 0.5, 0.5,
                           0.5, sd_within^2, 0.5,
                           0.5, 0.5, sd_within^2),
                         ncol = 3, nrow = 3)

  if(shared_re){

    random_effects <- mvrnorm(J, c(0,0,0),
                              Sigma = Sigma_random)


    db[,4,] <- matrix(rep(random_effects[db[,1,1], 1], m), ncol = m)
    db[,5,] <- matrix(rep(random_effects[db[,1,1], 2], m), ncol = m)
    db[,6,] <- matrix(rep(random_effects[db[,1,1], 3], m), ncol = m)


  }else{

    for(i in seq_len(m)){
      random_effects <- mvrnorm(J, c(0,0,0),
                                Sigma = Sigma_random)

      db[,4,i] <- random_effects[db[,1,1], 1]
      db[,5,i] <- random_effects[db[,1,1], 2]
      db[,6,i] <- random_effects[db[,1,1], 3]

    }


  }

  db_list <- list()
  for( i in seq_len(m)){

    db_list[[i]] <- data.frame(db[,,i])
    colnames(db_list[[i]]) <- c("id", "X", "Z", "U", "RX", "RZ", "eps")
    db_list[[i]]$y <- db_list[[i]]$X * beta[i] + db_list[[i]]$Z * gamma+
      db_list[[i]]$U + db_list[[i]]$RX * db_list[[i]]$X +
      db_list[[i]]$RZ * db_list[[i]]$Z  + db_list[[i]]$eps
  }

  attr(db_list, "Sigma_eps_true") <- Sigma
  attr(db_list, "Sigma_random_true")   <- Sigma_random

  return(db_list)
}
