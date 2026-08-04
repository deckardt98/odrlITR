test_that("odrl exposes parametric nuisances and user outer folds", {
  set.seed(84)
  n <- 90
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  a <- rep(c(-1, 1), length.out = n)
  y <- 1 + x[, 1] + a * x[, 2] + rnorm(n, sd = 0.2)
  fold_id <- rep(c("site-a", "site-b", "site-c"), length.out = n)

  fit <- odrl(
    x, a, y,
    learner = "relu", loss = "squared_hinge",
    nuisance = "parametric", nuisance_fold_id = fold_id,
    control = odrl_control(
      relu_architectures = list(integer()),
      relu_activation = "linear", relu_decay = 0.01,
      relu_folds = 2, relu_restarts = 1,
      relu_refit_restarts = 1, relu_maxit = 150,
      seed = 10
    )
  )

  expect_s3_class(fit, "odrl_fit")
  expect_identical(fit$nuisance$engine, "parametric")
  expect_identical(fit$nuisance$fold_source, "user supplied")
  expect_identical(fit$nuisance$fold_labels,
                   c("site-a", "site-b", "site-c"))
  expect_identical(fit$loss, "squared_hinge")
})

test_that("odrl exposes flexible SVM kernels and losses", {
  set.seed(85)
  n <- 48
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  a <- rep(c(-1, 1), length.out = n)
  y <- a * (x[, 1] - 0.5 * x[, 2]) + rnorm(n, sd = 0.1)
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = rep(0.5, n), out_of_fold = TRUE
  )

  fit <- odrl(
    x, a, y,
    learner = "svm", loss = "squared_hinge", nuisance = nuisance,
    control = odrl_control(
      svm_kernel = "linear", svm_hinge_mode = "regularized",
      svm_penalty = 0.1, svm_folds = 2, svm_maxit = 200,
      seed = 11
    )
  )

  expect_s3_class(fit, "odrl_fit")
  expect_identical(fit$policy$kernel, "linear")
  expect_identical(fit$loss, "squared_hinge")
  expect_length(predict(fit, x), n)
})

test_that("odrl adapts custom losses for neural policies", {
  set.seed(86)
  n <- 48
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  a <- rep(c(-1, 1), length.out = n)
  y <- a * (x[, 1] + 0.25) + rnorm(n, sd = 0.1)
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = rep(0.5, n), out_of_fold = TRUE
  )
  quadratic <- list(
    name = "quadratic_margin",
    value = function(margin) (1 - margin)^2,
    gradient = function(margin) -2 * (1 - margin)
  )

  fit <- odrl(
    x, a, y,
    learner = "relu", loss = quadratic, nuisance = nuisance,
    control = odrl_control(
      relu_architectures = list(integer()),
      relu_activation = "linear", relu_decay = 0.01,
      relu_folds = 2, relu_restarts = 1,
      relu_refit_restarts = 1, relu_maxit = 150,
      seed = 12
    )
  )

  expect_s3_class(fit, "odrl_fit")
  expect_identical(fit$loss, "quadratic_margin")
  expect_identical(fit$policy$diagnostics$architecture, integer())
  expect_identical(fit$policy$diagnostics$activation, "linear")
})

test_that("sparse scores do not break neural cross-validation", {
  x <- matrix(seq(-1, 1, length.out = 12), ncol = 2)
  a <- rep(c(1, -1), 3)
  y <- c(1, rep(0, 5))
  nuisance <- odrl_nuisance_user(
    m = rep(0, 6), pi = rep(0.5, 6), out_of_fold = TRUE
  )
  fit <- odrl(
    x, a, y, learner = "relu", loss = "logistic", nuisance = nuisance,
    control = odrl_control(
      relu_architectures = list(integer()), relu_activation = "linear",
      relu_decay = 0.01, relu_folds = 3, relu_restarts = 1,
      relu_refit_restarts = 1, relu_maxit = 50, seed = 13
    )
  )
  expect_s3_class(fit, "odrl_fit")
  expect_length(predict(fit, x), 6)
})
