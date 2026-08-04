test_that("built-in exponential loss is exact on ordinary margins", {
  margin <- c(-3, -1, 0, 0.5, 2, 5)
  specification <- odrlITR:::.odrl_resolve_surrogate_loss("exponential")

  expect_true(specification$builtin)
  expect_identical(specification$name, "exponential")
  expect_equal(specification$value(margin), exp(-margin), tolerance = 1e-14)
  expect_equal(
    specification$gradient(margin), -exp(-margin), tolerance = 1e-14
  )
})

test_that("exponential continuation is finite, smooth, and convex", {
  cutoff <- odrlITR:::.odrl_exponential_cutoff
  evaluate <- odrlITR:::.odrl_exponential_margin

  at_cutoff <- evaluate(-cutoff)
  expect_equal(at_cutoff$loss, exp(cutoff))
  expect_equal(at_cutoff$gradient, -exp(cutoff))

  epsilon <- 1e-7
  left <- evaluate(-cutoff - epsilon)
  right <- evaluate(-cutoff + epsilon)
  expect_equal(left$loss / at_cutoff$loss, 1 + epsilon, tolerance = 1e-12)
  expect_equal(
    right$loss / at_cutoff$loss, exp(-epsilon), tolerance = 1e-12
  )
  expect_equal(left$gradient / at_cutoff$gradient, 1, tolerance = 1e-12)
  expect_equal(
    right$gradient / at_cutoff$gradient, exp(-epsilon), tolerance = 1e-12
  )

  tail_margin <- -cutoff - 5
  step <- 1e-5
  numeric_gradient <- (
    evaluate(tail_margin + step)$loss -
      evaluate(tail_margin - step)$loss
  ) / (2 * step)
  expect_equal(
    numeric_gradient / evaluate(tail_margin)$gradient,
    1,
    tolerance = 1e-6
  )

  extreme <- evaluate(c(-1e6, 1e6))
  expect_true(all(is.finite(extreme$loss)))
  expect_true(all(is.finite(extreme$gradient)))
  expect_true(all(extreme$loss >= 0))
  expect_true(all(diff(evaluate(c(-40, -35, -30, -20, 0, 20))$gradient) >= 0))
})

test_that("SVM and native neural learners fit exponential loss", {
  set.seed(91)
  x <- matrix(rnorm(80), ncol = 2)
  score <- x[, 1] - 0.4 * x[, 2] + 0.1

  svm_control <- odrl_control(
    svm_kernel = "linear",
    svm_penalty = 0.1,
    svm_folds = 2,
    svm_maxit = 200,
    seed = 15
  )
  svm <- odrlITR:::.odrl_fit_svm(x, score, svm_control, "exponential")
  expect_identical(svm$loss, "exponential")
  expect_identical(svm$kernel, "linear")
  expect_false(svm$diagnostics$custom_loss)
  expect_false(svm$diagnostics$globally_bounded)
  expect_true("lambda" %in% names(svm$selected))
  expect_identical(svm$fit$convergence, 0L)
  expect_true(is.finite(svm$fit$objective))
  expect_true(all(is.finite(svm$tuning$mean_criterion)))
  expect_true(is.finite(svm$selected$mean_criterion))
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_svm(svm, x, type = "score")
  )))

  neural <- odrlITR:::.odrl_fit_relu_candidate(
    x, score, hidden = integer(), decay = 0.01,
    loss = "exponential", restarts = 1, maxit = 200, seed = 16,
    activation = "linear"
  )
  expect_identical(neural$loss_name, "exponential")
  expect_false(neural$bounded_output)
  expect_identical(neural$convergence, 0L)
  expect_true(neural$loss_builtin)
  expect_true(is.finite(neural$objective))
  expect_true(all(is.finite(neural$all_objectives)))
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_relu_raw(neural, x)
  )))

  prepared <- odrlITR:::.odrl_prepare_loss("relu", "exponential")
  expect_identical(prepared$fit, "exponential")
  expect_identical(prepared$label, "exponential")
})

test_that("native neural exponential gradient matches finite differences", {
  set.seed(92)
  x <- matrix(rnorm(18), nrow = 6, ncol = 3)
  label <- rep(c(-1, 1), 3)
  weight <- seq(0.6, 1.4, length.out = 6)
  architecture <- c(3L, 2L)
  theta <- odrlITR:::.odrl_relu_initial(
    x, label, architecture, seed = 17, loss = "exponential",
    activation = "tanh", weight = weight
  )
  objective <- function(parameter) {
    odrlITR:::.odrl_relu_objective(
      parameter, x, label, weight, architecture, decay = 0.02,
      loss = "exponential", activation = "tanh"
    )
  }
  analytic <- odrlITR:::.odrl_relu_objective(
    theta, x, label, weight, architecture, decay = 0.02,
    loss = "exponential", gradient = TRUE, activation = "tanh"
  )
  step <- 1e-6
  numeric_gradient <- vapply(seq_along(theta), function(index) {
    plus <- minus <- theta
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (objective(plus) - objective(minus)) / (2 * step)
  }, numeric(1))
  expect_equal(analytic, numeric_gradient, tolerance = 2e-5)
})

test_that("built-in and custom exponential losses agree on ordinary margins", {
  margin <- c(-4, -1, 0, 1, 4)
  builtin <- odrlITR:::.odrl_resolve_surrogate_loss("exponential")
  custom <- odrlITR:::.odrl_resolve_surrogate_loss(list(
    name = "custom_exponential",
    value = function(value) exp(-value),
    gradient = function(value) -exp(-value)
  ))
  expect_equal(builtin$value(margin), custom$value(margin))
  expect_equal(builtin$gradient(margin), custom$gradient(margin))

  set.seed(93)
  x <- matrix(rnorm(12), nrow = 4, ncol = 3)
  label <- c(-1, 1, -1, 1)
  weight <- c(0.7, 1.2, 0.9, 1.1)
  theta <- odrlITR:::.odrl_relu_initial(
    x, label, integer(), seed = 18, loss = "exponential",
    activation = "linear", weight = weight
  )
  callback <- function(value) {
    list(loss = exp(-value), gradient = -exp(-value))
  }
  expect_equal(
    odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, integer(), 0.01, "exponential"
    ),
    odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, integer(), 0.01, callback
    )
  )
  expect_equal(
    odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, integer(), 0.01, "exponential",
      gradient = TRUE
    ),
    odrlITR:::.odrl_relu_objective(
      theta, x, label, weight, integer(), 0.01, callback,
      gradient = TRUE
    )
  )
})

test_that("custom losses named exponential remain custom", {
  custom <- list(
    name = "exponential",
    value = function(margin) (1 - margin)^2,
    gradient = function(margin) -2 * (1 - margin)
  )
  specification <- odrlITR:::.odrl_resolve_surrogate_loss(custom)
  expect_false(specification$builtin)
  expect_equal(specification$value(c(-1, 0, 1)), c(4, 1, 0))
  expect_equal(specification$gradient(c(-1, 0, 1)), c(-4, -2, 0))

  neural <- odrlITR:::.odrl_prepare_loss("relu", custom)
  expect_true(is.function(neural$fit))
  expect_identical(neural$label, "exponential")
  evaluated <- neural$fit(c(-1, 0, 1))
  expect_equal(evaluated$loss, c(4, 1, 0))
  expect_equal(evaluated$gradient, c(-4, -2, 0))
})
