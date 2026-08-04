test_that("neural presets are transparent and explicit controls win", {
  standard <- odrl_control(relu_preset = "standard")
  expect_identical(
    standard$relu_architectures,
    list(integer(), 8L, 16L)
  )
  expect_identical(standard$relu_preset, "standard")
  expect_identical(standard$relu_backend, "native")

  overridden <- odrl_control(
    relu_preset = "flexible",
    relu_architectures = list(3L),
    relu_activation = "sigmoid",
    relu_restarts = 1L
  )
  expect_identical(overridden$relu_architectures, list(3L))
  expect_identical(overridden$relu_activation, "sigmoid")
  expect_identical(overridden$relu_restarts, 1L)
  expect_identical(overridden$relu_backend, "native")

  nnet <- odrl_control(relu_preset = "nnet")
  expect_identical(nnet$relu_backend, "nnet")
  expect_identical(nnet$relu_activation, "sigmoid")
  expect_identical(nnet$relu_architectures, list(integer(), 4L, 8L, 16L))
  expect_true(nnet$relu_backend_options$skip)

  direct_nnet <- odrl_control(relu_backend = "nnet")
  expect_identical(direct_nnet$relu_activation, "sigmoid")
  expect_error(
    odrl_control(relu_backend = "nnet", relu_activation = "relu"),
    "requires.*sigmoid"
  )
  expect_error(
    odrl_control(
      relu_backend = "nnet", relu_architectures = list(c(8L, 4L))
    ),
    "at most one hidden layer"
  )
})

test_that("finite-series specifications run through the public SVM route", {
  set.seed(931)
  n <- 60
  x <- cbind(x1 = runif(n), x2 = runif(n))
  a <- rep(c(-1, 1), length.out = n)
  y <- a * (x[, 1] - 0.4) + rnorm(n, sd = 0.1)
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = rep(0.5, n), out_of_fold = TRUE
  )
  specification <- odrl_series_kernel(
    "legendre", legendre_degree = c(1L, 2L),
    combine = "additive", bounds = c(0, 1)
  )
  fit <- odrl(
    x, a, y, learner = "svm", loss = "logistic", nuisance = nuisance,
    control = odrl_control(
      svm_kernel = specification, svm_penalty = 0.1,
      svm_folds = 2L, svm_maxit = 300L, seed = 932L
    )
  )

  expect_identical(fit$policy$engine, "series")
  expect_identical(fit$policy$kernel, "legendre")
  expect_true(fit$policy$diagnostics$finite_series)
  expect_gt(fit$policy$diagnostics$series_features, 0)
  expect_null(fit$policy$transform)
  expect_length(predict(fit, x), n)
})

test_that("the nnet quick-start backend fits weighted logistic policies", {
  skip_if_not_installed("nnet")
  set.seed(933)
  n <- 60
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  a <- rep(c(-1, 1), length.out = n)
  y <- a * (x[, 1] - 0.25 * x[, 2]) + rnorm(n, sd = 0.1)
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = rep(0.5, n), out_of_fold = TRUE
  )
  fit <- odrl(
    x, a, y, learner = "relu", loss = "logistic", nuisance = nuisance,
    control = odrl_control(
      relu_preset = "nnet", relu_architectures = list(integer()),
      relu_decay = 0.001, relu_folds = 2L, relu_restarts = 1L,
      relu_refit_restarts = 1L, relu_maxit = 200L, seed = 934L
    )
  )

  expect_identical(fit$policy$backend, "nnet")
  expect_identical(fit$policy$diagnostics$preset, "nnet")
  expect_identical(fit$policy$diagnostics$activation, "sigmoid")
  expect_length(predict(fit, x), n)
})

test_that("automatic hinge mode preserves the special RBF route only", {
  rbf <- odrl_control(svm_kernel = "rbf")
  linear <- odrl_control(svm_kernel = "linear")
  expect_identical(
    odrlITR:::.odrl_effective_hinge_mode(
      odrlITR:::.odrl_resolve_kernel(rbf$svm_kernel, rbf),
      rbf$svm_hinge_mode
    ),
    "bounded"
  )
  expect_identical(
    odrlITR:::.odrl_effective_hinge_mode(
      odrlITR:::.odrl_resolve_kernel(linear$svm_kernel, linear),
      linear$svm_hinge_mode
    ),
    "regularized"
  )
})
