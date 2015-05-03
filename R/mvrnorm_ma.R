#' Simulate from a Multivariate Normal Distribution with MA(q) dependence
#'
#' @description
#' Produces one or more samples from the specified multivariate normal distribution,
#' in which the errors are MA(q).
#'
#' @param n the number of samples required.
#' @param mu	a vector giving the means of the variables.
#' @param Sigma	a positive-definite symmetric matrix specifying the covariance matrix of the variables.
#' @param tol	tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
#' @param empirical	logical. If true, mu and Sigma specify the empirical not population mean and covariance matrix.
#' @param rho A numeric vector of MA coeficients.
#' @return A matrix
#' @export
mvrnorm_ma <- function (n = 1, mu, Sigma, rho = 0, tol = 1e-06, empirical = FALSE)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(0, nrow = n, ncol = p)
  Xerr <- rep(0, p)
  for (i in seq_len(n)) {
    Xerr <- rnorm(p) + rho * Xerr
    X[i, ] <- Xerr
  }
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}
