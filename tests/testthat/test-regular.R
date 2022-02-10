set.seed(42)
library(spectralGraphTopology)

test_that("test learn_connected_graph", {
  w <- c(1, 1, 1, 1, 1, 1) / 3
  Laplacian <- L(w)
  p <- ncol(Laplacian)
  S <- cov(MASS::mvrnorm(p * 100, rep(0, p), MASS::ginv(Laplacian)))
  res <- learn_connected_graph(S, rho = 100)
  laplacian <- res$laplacian
  expect_true(res$convergence)
  expect_true(res$maxiter > 5)
  expect_true(spectralGraphTopology:::relative_error(Laplacian, laplacian) < 1e-1)
  expect_true(spectralGraphTopology:::fscore(Laplacian, laplacian, 1e-1) > .9)
})

test_that("test learn_regular_heavytail_graph", {
  w <- c(1, 1, 1, 1, 1, 1) / 3
  Laplacian <- L(w)
  p <- ncol(Laplacian)
  X <- MASS::mvrnorm(p * 100, rep(0, p), MASS::ginv(Laplacian))
  res <- learn_regular_heavytail_graph(X, rho = 100, heavy_type = "student", nu = 1e3)
  laplacian <- res$laplacian
  expect_true(res$convergence)
  expect_true(res$maxiter > 5)
  expect_true(relative_error(Laplacian, laplacian) < 1e-1)
  expect_true(fscore(Laplacian, laplacian, 1e-1) > .9)
})

test_that("test learn_regular_heavytail_graph", {
  w <- c(1, 1, 1, 1, 1, 1) / 3
  Laplacian <- L(w)
  p <- ncol(Laplacian)
  X <- scale(MASS::mvrnorm(p * 100, rep(0, p), MASS::ginv(Laplacian)))
  res_1 <- learn_regular_heavytail_graph(X, rho = 100)
  res_2 <- learn_regular_heavytail_graph(X, rho = 100, heavy_type = "student", nu = 1e4)
  res_3 <- learn_connected_graph(cor(X), rho = 100)
  expect_true(relative_error(res_1$laplacian, res_2$laplacian) < 1e-4)
  expect_true(relative_error(res_2$laplacian, res_3$laplacian) < 1e-4)
})

test_that("test learn_regular_heavytail_graph student", {
  w <- c(1, 1, 1, 1, 1, 1) / 3
  Laplacian <- L(w)
  p <- ncol(Laplacian)
  nu <- 4
  X <- mvtnorm::rmvt(n = p * 500, delta = rep(0, p), sigma = ((nu-2)/nu) * MASS::ginv(Laplacian), df = nu)
  res <- learn_regular_heavytail_graph(X, rho = 1, heavy_type = "student", nu = nu)
  laplacian <- res$laplacian
  expect_true(res$convergence)
  expect_true(res$maxiter > 5)
  expect_true(relative_error(Laplacian, laplacian) < 1e-1)
  expect_true(fscore(Laplacian, laplacian, 1e-2) > .9)
})


test_that("test learn_kcomp_heavytail_graph", {
  w1 <- c(1, 1, 1, 1, 1, 1)/3
  w2 <- c(1, 1, 1, 1, 1, 1)/3
  Laplacian <- block_diag(L(w1), L(w2))
  p <- ncol(Laplacian)
  nu <- 4
  X <- mvtnorm::rmvt(n = p * 500, delta = rep(0, p), sigma = ((nu-2)/nu) * MASS::ginv(Laplacian), df = nu)
  res <- learn_kcomp_heavytail_graph(X, k = 2, rho = 1e2, heavy_type = "student", nu = nu, reltol = 1e-4)
  laplacian <- res$laplacian
  expect_true(res$convergence)
  expect_true(res$maxiter > 5)
  expect_true(spectralGraphTopology:::relative_error(Laplacian, laplacian) < 1e-1)
  expect_true(spectralGraphTopology:::fscore(Laplacian, laplacian, 1e-2) > .9)
})


test_that("test learn_kcomp_heavytail_graph", {
  w1 <- c(1, 1, 1, 1, 1, 1)/3
  w2 <- c(1, 1, 1, 1, 1, 1)/3
  Laplacian <- block_diag(L(w1), L(w2))
  p <- ncol(Laplacian)
  X <- MASS::mvrnorm(p * 100, rep(0, p), MASS::ginv(Laplacian))
  res <- learn_kcomp_heavytail_graph(X, k = 2, rho = 100)
  laplacian <- res$laplacian
  expect_true(res$convergence)
  expect_true(res$maxiter > 5)
  expect_true(spectralGraphTopology:::relative_error(Laplacian, laplacian) < 1e-1)
  expect_true(spectralGraphTopology:::fscore(Laplacian, laplacian, 1e-1) > .9)
})
