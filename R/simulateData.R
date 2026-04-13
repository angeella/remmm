#' @title Simulate multivariate hierarchical data
#'
#' @description
#' Generates multivariate hierarchical data with \eqn{m} outcomes, shared
#' cluster-level random effects, and correlated residual errors across outcomes.
#'
#' For each outcome \eqn{i = 1, \dots, m}, the response is generated as:
#'
#' \deqn{
#' y_i = X \beta_i + Z \gamma + U + RX \cdot X + RZ \cdot Z + \varepsilon_i
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{X} and \eqn{Z} are shared fixed-effect covariates,
#'   \item \eqn{U}, \eqn{RX}, \eqn{RZ} are shared cluster-level random effects,
#'   \item \eqn{\varepsilon} follows a multivariate normal distribution
#'         with covariance matrix \code{Sigma}.
#' }
#'
#' @param Sigma Either:
#'   \itemize{
#'     \item \code{"autocorrelation"}: AR(1)-like structure with correlation \code{rho},
#'     \item \code{"equicorrelation"}: constant correlation \code{rho},
#'     \item or a numeric covariance matrix.
#'   }
#'   Default \code{"equicorrelation"}.
#' @param rho Correlation parameter used when \code{Sigma} is
#'   \code{"autocorrelation"} or \code{"equicorrelation"}. Default \code{0.5}.
#' @param beta Numeric vector of length \eqn{m}; outcome-specific slopes for \code{X}.
#' @param gamma Numeric scalar; common slope for \code{Z} across outcomes.
#' @param J Integer; number of clusters.
#' @param nJ Integer vector of length \code{J}; number of observations per cluster.
#' @param sd_within Numeric; standard deviation of cluster-level random effects.
#'   Default \code{1}.
#' @param seed Integer; random seed for reproducibility. Default \code{1234}.
#' @param sd_eps Numeric; marginal standard deviation of the residual errors.
#'   Default \code{1}.
#'
#' @return A single \code{data.frame} with columns:
#' \describe{
#'   \item{\code{id}}{Cluster identifier}
#'   \item{\code{X}, \code{Z}}{Shared fixed-effect covariates}
#'   \item{\code{U}, \code{RX}, \code{RZ}}{Shared random effects}
#'   \item{\code{eps}}{Matrix column of residual terms, one column per outcome}
#'   \item{\code{Y}}{Matrix column of responses, one column per outcome}
#' }
#'
#' The returned data frame contains two attributes:
#' \describe{
#'   \item{\code{Sigma_eps_true}}{True residual covariance matrix}
#'   \item{\code{Sigma_random_true}}{True random-effects covariance matrix}
#' }
#'
#' @examples
#' db <- simulateData(Sigma = "equicorrelation", rho = 0.5, beta = c(0.5, 1.2), gamma = 0.8,
#'   J = 8, nJ = rep(25, 8))
#'
#' str(db)
#'
#' @importFrom MASS mvrnorm
#' @export

simulateData <- function(Sigma = "equicorrelation", rho = 0.5,
                         beta, gamma,
                         J, nJ, sd_within = 1,
                         seed = 1234, sd_eps = 1) {

  set.seed(seed)

  m <- length(beta)
  N <- sum(nJ)

  if (length(nJ) != J) {
    stop("length(nJ) must be equal to J")
  }

  if (!is.matrix(Sigma)) {

    if (!Sigma %in% c("autocorrelation", "equicorrelation")) {
      stop("Sigma must be either a covariance matrix, 'autocorrelation', or 'equicorrelation'")
    }

    if (Sigma == "autocorrelation") {
      if (abs(rho) >= 1) {
        stop("For 'autocorrelation', rho must satisfy |rho| < 1")
      }
      Sigma <- sd_eps^2 * outer(1:m, 1:m, function(i, j) rho^abs(i - j))
    }

    if (Sigma == "equicorrelation") {
      if (rho <= -1 / (m - 1) || rho >= 1) {
        stop("For 'equicorrelation', rho must satisfy -1/(m-1) < rho < 1")
      }
      Sigma <- matrix(rho * sd_eps^2, m, m)
      diag(Sigma) <- sd_eps^2
    }

  } else {
    if (!all(dim(Sigma) == c(m, m))) {
      stop("Numeric Sigma must be an m x m covariance matrix")
    }
  }

  id <- rep(seq_len(J), times = nJ)

  Sigma_fixed <- matrix(c(1, 0.7,
                          0.7, 1), ncol = 2)

  fixed_effects <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_fixed)
  X <- fixed_effects[, 1]
  Z <- fixed_effects[, 2]

  E <- mvrnorm(N, mu = rep(0, m), Sigma = Sigma)
  colnames(E) <- paste0("eps", seq_len(m))

  R_random <- matrix(0.5, 3, 3)
  diag(R_random) <- 1
  Sigma_random <- sd_within^2 * R_random

  random_effects <- mvrnorm(J, mu = c(0, 0, 0), Sigma = Sigma_random)

  U  <- random_effects[id, 1]
  RX <- random_effects[id, 2]
  RZ <- random_effects[id, 3]

  mean_common <- U + RX * X + RZ * Z + gamma * Z
  Y <- sweep(X %o% beta + E, 1, mean_common, FUN = "+")
  colnames(Y) <- paste0("y", seq_len(m))

  db <- data.frame(
    id = id,
    X = X,
    Z = Z,
    U = U,
    RX = RX,
    RZ = RZ
  )

  db$eps <- I(E)
  db$Y <- I(Y)

  attr(db, "Sigma_eps_true") <- Sigma
  attr(db, "Sigma_random_true") <- Sigma_random

  db
}
