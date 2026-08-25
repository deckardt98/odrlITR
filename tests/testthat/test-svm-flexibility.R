.svm_flex_control <- function(kernel = "rbf") {
  control <- odrl_control(
    svm_kernel = if (is.character(kernel) &&
      kernel %in% c("rbf", "linear")) kernel else "rbf",
    svm_rbf_multiplier = 1,
    svm_penalty = 0.1,
    svm_folds = 2,
    svm_maxit = 200,
    seed = 31
  )
  control$svm_kernel <- kernel
  control$svm_hinge_mode <- "bounded"
  control$svm_polynomial_degree <- 2L
  control$svm_polynomial_scale <- 1
  control$svm_polynomial_offset <- 1
  control
}

test_that("common and custom kernel matrices are available", {
  x <- matrix(c(-1, 0, 1, 1, 0, -1), ncol = 2)
  linear <- odrlITR:::.odrl_kernel_matrix(x, kernel = "linear")
  gaussian <- odrlITR:::.odrl_kernel_matrix(
    x, kernel = "gaussian", bandwidth2 = 1
  )
  polynomial <- odrlITR:::.odrl_kernel_matrix(
    x,
    kernel = list(
      name = "polynomial", degree = 2L, scale = 1, offset = 1
    )
  )
  laplacian <- function(x, y) {
    distance <- sqrt(odrlITR:::.odrl_squared_distance(x, y))
    exp(-distance)
  }
  custom <- odrlITR:::.odrl_kernel_matrix(x, kernel = laplacian)
  custom_named_like_builtin <- odrlITR:::.odrl_kernel_matrix(
    x,
    kernel = list(
      name = "linear",
      fun = function(x, y) matrix(0.25, nrow(x), nrow(y))
    )
  )

  expect_equal(dim(linear), c(3L, 3L))
  expect_equal(diag(gaussian), rep(1, 3))
  expect_equal(polynomial, (linear + 1)^2)
  expect_equal(diag(custom), rep(1, 3))
  expect_equal(custom, t(custom), tolerance = 1e-12)
  expect_equal(custom_named_like_builtin, matrix(0.25, 3, 3))
})

test_that("squared hinge works with linear and polynomial kernels", {
  set.seed(41)
  x <- matrix(rnorm(80), ncol = 2)
  score <- x[, 1] - 0.4 * x[, 2] + 0.05

  linear_control <- .svm_flex_control("linear")
  linear_control$svm_hinge_mode <- "regularized"
  linear <- odrlITR:::.odrl_fit_svm(
    x, score, linear_control, "squared_hinge"
  )

  polynomial_control <- .svm_flex_control("polynomial")
  polynomial_control$svm_polynomial_degree <- c(2L, 3L)
  polynomial <- odrlITR:::.odrl_fit_svm(
    x, score, polynomial_control, "squared_hinge"
  )

  expect_s3_class(linear, "odrl_policy_svm")
  expect_identical(linear$loss, "squared_hinge")
  expect_identical(linear$kernel, "linear")
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_svm(linear, x, type = "score")
  )))
  expect_identical(polynomial$kernel, "polynomial")
  expect_true(polynomial$selected$degree %in% c(2L, 3L))
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_svm(polynomial, x, type = "score")
  )))
})

test_that("bounded linear hinge fits first and clips second", {
  x <- sqrt(2) * diag(2)
  score <- c(1, -1)
  control <- odrl_control(
    svm_kernel = "linear", svm_hinge_mode = "bounded",
    svm_penalty = 1, svm_folds = 2, svm_maxit = 500
  )
  fit <- odrlITR:::.odrl_fit_svm_candidate(
    x, score, control, "hinge",
    candidate = list(multiplier = 1, lambda = 1), seed = 47
  )

  expect_equal(fit$intercept, 0, tolerance = 1e-6)
  expect_equal(fit$alpha, c(0.5, -0.5), tolerance = 1e-5)
  expect_equal(fit$objective, 0.75, tolerance = 1e-6)
  expect_equal(fit$rkhs_norm, sqrt(0.5), tolerance = 1e-5)
  expect_identical(fit$clipping, "hard_tanh")
  expect_true(fit$globally_bounded)

  newx <- matrix(c(10 * sqrt(2), 0), nrow = 1)
  raw <- odrlITR:::.odrl_predict_kernel_unclipped(fit, newx)
  bounded <- odrlITR:::.odrl_predict_kernel_raw(fit, newx)
  expect_equal(raw, 5, tolerance = 1e-5)
  expect_equal(bounded, 1)
  expect_equal(sign(raw), sign(bounded))

  kink_fit <- odrlITR:::.odrl_fit_svm_candidate(
    x, score, control, "hinge",
    candidate = list(multiplier = 1, lambda = 0.01), seed = 48
  )
  expect_equal(kink_fit$convergence, 0)
  expect_equal(kink_fit$alpha, c(1, -1), tolerance = 1e-5)
  expect_equal(kink_fit$objective, 0.01, tolerance = 1e-6)
  expect_lt(kink_fit$duality_gap, 1e-6)
})

test_that("custom kernel and signed-margin loss are validated and retained", {
  set.seed(52)
  x <- matrix(rnorm(72), ncol = 2)
  score <- x[, 1] * x[, 2] + 0.1
  laplacian <- list(
    name = "laplacian",
    fun = function(x, y, rate) {
      distance <- sqrt(odrlITR:::.odrl_squared_distance(x, y))
      exp(-rate * distance)
    },
    args = list(rate = 0.75)
  )
  smooth_margin <- list(
    name = "softplus_margin",
    value = function(margin) odrlITR:::.odrl_log1pexp(1 - margin),
    gradient = function(margin) -stats::plogis(1 - margin)
  )
  control <- .svm_flex_control(laplacian)
  policy <- odrlITR:::.odrl_fit_svm(x, score, control, smooth_margin)

  expect_identical(policy$kernel, "laplacian")
  expect_identical(policy$loss, "softplus_margin")
  expect_identical(policy$diagnostics$kernel, "laplacian")
  expect_identical(policy$diagnostics$loss, "softplus_margin")
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_svm(policy, x, type = "score")
  )))
})

test_that("custom loss labels never select built-in hinge semantics", {
  set.seed(53)
  x <- matrix(rnorm(60), ncol = 2)
  score <- x[, 1] - x[, 2]
  custom_hinge_label <- list(
    name = "hinge",
    value = function(margin) (1 - margin)^2,
    gradient = function(margin) -2 * (1 - margin)
  )
  control <- .svm_flex_control("rbf")
  grid <- odrlITR:::.odrl_kernel_grid(control)
  expect_true("lambda" %in% names(grid))
  expect_false("radius" %in% names(grid))

  policy <- odrlITR:::.odrl_fit_svm(
    x, score, control, custom_hinge_label
  )
  expect_true(policy$diagnostics$custom_loss)
  expect_false(policy$diagnostics$globally_bounded)
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_svm(policy, x, type = "score")
  )))
})

test_that("invalid custom SVM specifications fail clearly", {
  x <- matrix(rnorm(20), ncol = 2)
  expect_error(
    odrlITR:::.odrl_kernel_matrix(
      x, kernel = function(x, y) matrix(1, 1, 1)
    ),
    "nrow\\(x\\).*nrow\\(y\\)"
  )
  expect_error(
    odrlITR:::.odrl_resolve_surrogate_loss(
      list(value = function(margin) margin^2)
    ),
    "gradient"
  )
  expect_error(
    odrlITR:::.odrl_resolve_surrogate_loss(
      list(
        value = function(margin) rep(-1, length(margin)),
        gradient = function(margin) rep(0, length(margin))
      )
    ),
    "nonnegative"
  )
})
