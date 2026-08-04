test_that("native neural backend supports deep architectures and activations", {
  set.seed(41)
  x <- matrix(rnorm(24), nrow = 8, ncol = 3)
  score <- rep(c(-1, 1), 4)
  label <- sign(score)
  weight <- rep(1, length(score))
  architecture <- c(4L, 2L)

  for (activation in c("relu", "leaky_relu", "tanh", "sigmoid", "linear")) {
    theta <- odrlITR:::.odrl_relu_initial(
      x, label, architecture, seed = 12, loss = "squared_hinge",
      activation = activation
    )
    value <- odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, architecture, decay = 0.01,
      loss = "squared_hinge", activation = activation
    )
    gradient <- odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, architecture, decay = 0.01,
      loss = "squared_hinge", gradient = TRUE, activation = activation
    )
    expect_true(is.finite(value), info = activation)
    expect_equal(length(gradient), length(theta), info = activation)
    expect_true(all(is.finite(gradient)), info = activation)
  }
})

test_that("deep-network analytic gradient agrees with finite differences", {
  set.seed(42)
  x <- matrix(rnorm(18), nrow = 6, ncol = 3)
  label <- rep(c(-1, 1), 3)
  weight <- seq(0.7, 1.2, length.out = 6)
  architecture <- c(3L, 2L)
  theta <- odrlITR:::.odrl_relu_initial(
    x, label, architecture, seed = 13, loss = "logistic",
    activation = "tanh"
  )
  objective <- function(parameter) {
    odrlITR:::.odrl_relu_objective(
      parameter, x, label, weight, architecture, decay = 0.02,
      loss = "logistic", activation = "tanh"
    )
  }
  analytic <- odrlITR:::.odrl_relu_objective(
    theta, x, label, weight, architecture, decay = 0.02,
    loss = "logistic", gradient = TRUE, activation = "tanh"
  )
  epsilon <- 1e-6
  numeric_gradient <- vapply(seq_along(theta), function(index) {
    plus <- minus <- theta
    plus[[index]] <- plus[[index]] + epsilon
    minus[[index]] <- minus[[index]] - epsilon
    (objective(plus) - objective(minus)) / (2 * epsilon)
  }, numeric(1))
  expect_equal(analytic, numeric_gradient, tolerance = 2e-5)
})

test_that("custom differentiable neural losses are validated and fitted", {
  custom_loss <- function(margin) {
    list(loss = (1 - margin)^2, gradient = -2 * (1 - margin))
  }
  attr(custom_loss, "name") <- "quadratic_margin"
  x <- matrix(seq(-1, 1, length.out = 20), ncol = 1)
  score <- ifelse(x[, 1] >= 0, 1, -1)
  fit <- odrlITR:::.odrl_fit_relu_candidate(
    x, score, hidden = integer(), decay = 0.01, loss = custom_loss,
    restarts = 1, maxit = 100, seed = 7, activation = "linear"
  )
  expect_identical(fit$loss_name, "quadratic_margin")
  expect_true(all(is.finite(odrlITR:::.odrl_predict_relu_raw(fit, x))))

  bad_loss <- function(margin) rep(0, length(margin))
  expect_error(
    odrlITR:::.odrl_relu_objective(
      c(0, 0), x, sign(score), rep(1, nrow(x)), integer(), 0,
      bad_loss
    ),
    "must return a list"
  )

  negative_loss <- function(margin) {
    list(loss = rep(-1, length(margin)), gradient = rep(0, length(margin)))
  }
  expect_error(
    odrlITR:::.odrl_relu_objective(
      c(0, 0), x, sign(score), rep(1, nrow(x)), integer(), 0,
      negative_loss
    ),
    "nonnegative"
  )
})

test_that("a custom neural backend is genuinely callable", {
  backend <- list(
    name = "mock_linear_backend",
    fit = function(x, score, architecture, activation, decay, loss,
                   restarts, maxit, seed, leaky_slope) {
      list(beta = drop(crossprod(x, score)))
    },
    predict = function(model, newx) drop(newx %*% model$beta)
  )
  x <- cbind(x1 = c(-2, -1, 1, 2), x2 = c(0, 1, 0, 1))
  score <- c(-2, -1, 1, 2)
  fit <- odrlITR:::.odrl_fit_relu_candidate(
    x, score, hidden = c(3L, 2L), decay = 0, loss = "logistic",
    restarts = 1, maxit = 2, seed = 1, backend = backend
  )
  expect_identical(fit$backend_name, "mock_linear_backend")
  expect_equal(fit$architecture, c(3L, 2L))
  expect_equal(
    odrlITR:::.odrl_predict_relu_raw(fit, x),
    drop(x %*% crossprod(x, score))
  )
})

test_that("custom loss labels do not clip neural predictions", {
  custom_hinge_label <- function(margin) {
    list(loss = (1 - margin)^2, gradient = -2 * (1 - margin))
  }
  attr(custom_hinge_label, "name") <- "hinge"
  backend <- list(
    name = "constant_two",
    fit = function(x, score, architecture, activation, decay, loss,
                   restarts, maxit, seed, leaky_slope) list(),
    predict = function(model, newx) rep(2, nrow(newx))
  )
  x <- matrix(c(-1, 1), ncol = 1)
  fit <- odrlITR:::.odrl_fit_relu_candidate(
    x, c(-1, 1), hidden = integer(), decay = 0,
    loss = custom_hinge_label, restarts = 1, maxit = 1, seed = 1,
    backend = backend
  )
  expect_false(fit$bounded_output)
  expect_equal(odrlITR:::.odrl_predict_relu_raw(fit, x), c(2, 2))
})

test_that("all-zero neural training scores use a deterministic tie fit", {
  x <- matrix(seq(-1, 1, length.out = 8), ncol = 1)
  fit <- odrlITR:::.odrl_fit_relu_candidate(
    x, rep(0, nrow(x)), hidden = c(2L), decay = 0.01,
    loss = "logistic", restarts = 1, maxit = 20, seed = 1
  )
  expect_identical(fit$backend_name, "constant-zero-score")
  expect_equal(odrlITR:::.odrl_predict_relu_raw(fit, x), rep(0, nrow(x)))
})

test_that("custom exact-tree backends receive depth and options", {
  backend <- list(
    name = "mock_tree",
    fit = function(x, score, rewards, depth, min_node_size, split_step,
                   options) {
      list(cut = options$cut, depth = depth, score = score)
    },
    predict = function(model, newx) ifelse(newx[, 1] >= model$cut, 1, -1)
  )
  control <- odrl_control(tree_depth = 3, tree_min_node_size = 2)
  control$tree_backend <- backend
  control$tree_options <- list(cut = 0)
  x <- matrix(c(-2, -1, 1, 2), ncol = 1)
  fit <- odrlITR:::.odrl_fit_tree(x, c(-2, -1, 1, 2), control)
  expect_identical(fit$engine, "mock_tree")
  expect_equal(fit$fit$depth, 3L)
  expect_equal(fit$training_action, c(-1, -1, 1, 1))
  expect_equal(
    odrlITR:::.odrl_predict_tree(fit, matrix(c(-3, 3), ncol = 1)),
    c(-1, 1)
  )
  expect_identical(fit$diagnostics$options$cut, 0)
})

test_that("affine solver persists advanced stopping controls", {
  skip_if_not_installed("highs")
  control <- odrl_control(linear_time_limit = 10, linear_relative_gap = 0.1)
  control$linear_absolute_gap <- 0.1
  control$linear_node_limit <- 1000L
  control$linear_objective_target <- NULL
  control$linear_threads <- 1L
  control$linear_log_to_console <- FALSE
  control$linear_solver_options <- list()
  x <- matrix(c(-2, -1, 1, 2), ncol = 1)
  fit <- odrlITR:::.odrl_fit_linear(x, c(-2, -1, 1, 2), control)
  expect_identical(fit$diagnostics$requested_absolute_gap, 0.1)
  expect_identical(fit$diagnostics$node_limit, 1000L)
  expect_identical(fit$diagnostics$threads, 1L)
  expect_true(is.logical(fit$diagnostics$absolute_gap_met))

  control$linear_solver_options <- list(time_limit = 1)
  expect_error(
    odrlITR:::.odrl_linear_solver_control(control, 1e-8),
    "cannot override"
  )
})
