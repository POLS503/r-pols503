#' @import tidyr
#' @import lintr
#' @import knitr
#' @import dplyr
#' @import lazyeval
NULL

# From tidyr:::append_df
append_df <- function (x, values, after = length(x)) {
  y <- append(x, values, after = after)
  class(y) <- class(x)
  attr(y, "row.names") <- attr(x, "row.names")
  y
}

# From tidyr:::append_col
append_col <- function (x, col, name, after = length(x)) {
  append_df(x, setNames(list(col), name), after = after)
}

# From tidyr:::col_name
col_name <- function (x, default = stop("Please supply column name",
                      call. = FALSE)) {
  if (is.character(x))
    return(x)
  if (identical(x, quote(expr = )))
    return(default)
  if (is.name(x))
    return(as.character(x))
  if (is.null(x))
    return(x)
  stop("Invalid column specification", call. = FALSE)
}

#' Create covariance matrix from standard deviations and correlation matrix
#'
#' A covariance matrix can be specified as,
#' \deqn{\Sigma = diag(s) R diag(s)}
#' where \eqn{\Sigma} is the covariance matrix, \eqn{diag(s)} is a diagonal matrix
#' with the standard deviation on its diagonal, and \eqn{R} is the correlation matrix.
#'
#' @param sd \code{numeric} A vector of standard deviations
#' @param cor \code{matrix} A correlation matrix.
#' @return A covariance matrix.
#' @export
#' @examples
#' s <- c(1, 2)
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
#' sdcor2cov(s, R)
sdcor2cov <- function(sd, cor = diag(1, nrow = length(sd), ncol = length(sd))) {
  sdiag <- diag(sd, nrow = length(sd), ncol = length(sd))
  sdiag %*% cov2cor(cor) %*% sdiag
}

cov2sdcor <- function(x) {
  list(cor = cov2cor(x), sd = sqrt(diag(x)))
}

#' Create an equicorrelated correlation matrix.
#'
#' A correlation matrix with all off diagonal correlations equal to \eqn{rho}.
#'
#' @param n Size of the matrix.
#' @param rho The off diagonal correlations
#' @return A correlation matrix.
#' @export
#' @examples
#' n <- 3
#' rho <- 0.25
#' equicor(n, rho)
equicor <- function(n, rho) {
  ret <- matrix(rho, nrow = n, ncol = n)
  diag(ret) <- rep(1, n)
  ret
}

#' Find regression standard error to produce a given R-squared
#'
#' @description For a normal linear model, given a sample X and population beta,
#' find the regression standard deviation to produce the desired R-squared.
#'
#' @param X The sample design matrix
#' @param b The population coefficients for X
#' @param r2 The desired R^2
#' @return The regression standard deviation to produce that R^2.
#' @export
r2_to_sigma <- function(X, b, r2) {
  n <- nrow(X)
  k <- ncol(X)
  yhat <- X %*% b
  # Var(Y) = E(Var(Y|X)) + Var(E(Y|X))
  ybar <- mean(yhat)
  # Model sum of squares
  ssm <- sum((yhat - ybar) ^ 2)
  # calc sum of squared errors
  sse <- (1 - r2) / r2 * ssm
  # R^2 = 1 - SSE / (SSM + SSE)
  # (1 - R^2) * (SSM + SSE) = SSE
  # (1 - R^2) * SSM = SSE - (1 - R^2) * SSE
  # ((1 - R^2) * SSM) / R^2 = SSE
  # return sigma (assume population n)
  sqrt(sse / n)
}
