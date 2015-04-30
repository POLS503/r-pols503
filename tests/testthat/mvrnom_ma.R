context("mvrnom_ma")

test_that("mvrnorm runs", {
  x <- mvrnorm_ma(100, mu = c(0, 0), Sigma = diag(2), rho = c(0.5, 0.2))
  expect_is(x, "matrix")
  expect_equivalent(dim(x), c(100, 2))
})

test_that("mvrnorm with empirical = TRUE returns exact equiv", {
  x <- mvrnorm_ma(100, mu = c(0, 1), Sigma = matrix(c(1, 0.5, 0.5, 1),
                                                    nrow = 2),
                  rho = 0.5, empirical = TRUE)
  expect_is(x, "matrix")
  expect_equivalent(dim(x), c(100, 2))
  expect_equivalent(cov(x), matrix(c(1, 0.5, 0.5, 1), nrow = 2))
  expect_equivalent(apply(x, 2, mean), c(0, 1))
})
