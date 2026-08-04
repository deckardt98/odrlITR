test_that("nnet backend options are validated", {
  options <- odrlITR:::.odrl_nnet_options()
  expect_true(options$skip)
  expect_false(options$trace)
  expect_equal(options$MaxNWts, 10000L)

  expect_error(
    odrlITR:::.odrl_nnet_options(list(unknown = 1)),
    "Unknown `nnet` backend option"
  )
  expect_error(
    odrlITR:::.odrl_nnet_options(list(skip = 1)),
    "must be TRUE or FALSE"
  )
  expect_error(
    odrlITR:::.odrl_nnet_options(list(probability_epsilon = 0.5)),
    "must lie"
  )
})

test_that("nnet backend rejects incompatible policy specifications", {
  skip_if_not_installed("nnet")
  x <- matrix(rnorm(40), ncol = 2)
  score <- rep(c(-1, 1), 10)
  common <- list(x = x, score = score, decay = 0, restarts = 1,
                 maxit = 20, seed = 1)

  expect_error(
    do.call(odrlITR:::.odrl_fit_nnet_backend,
            c(common, list(architecture = c(4L, 2L)))),
    "at most one hidden layer"
  )
  expect_error(
    do.call(odrlITR:::.odrl_fit_nnet_backend,
            c(common, list(architecture = 4.5))),
    "zero or one positive integer width"
  )
  expect_error(
    do.call(odrlITR:::.odrl_fit_nnet_backend,
            c(common, list(architecture = 4L, activation = "relu"))),
    "sigmoid hidden units"
  )
  for (loss in c("hinge", "exponential", "squared_hinge")) {
    expect_error(
      do.call(odrlITR:::.odrl_fit_nnet_backend,
              c(common, list(architecture = 4L, loss = loss))),
      "supports only.*logistic",
      info = loss
    )
  }
  expect_error(
    do.call(odrlITR:::.odrl_fit_nnet_backend,
            c(common, list(architecture = integer(),
                           options = list(skip = FALSE)))),
    "requires `skip = TRUE`"
  )
  expect_error(
    do.call(odrlITR:::.odrl_fit_nnet_backend,
            c(common, list(architecture = 20L,
                           options = list(MaxNWts = 10L)))),
    "needs.*weights"
  )
})

test_that("nnet backend fits affine and one-hidden-layer logistic scores", {
  skip_if_not_installed("nnet")
  set.seed(71)
  x <- matrix(rnorm(240), ncol = 2)
  score <- ifelse(x[, 1] - 0.4 * x[, 2] >= 0, 1, -1)

  affine <- odrlITR:::.odrl_fit_nnet_backend(
    x, score, architecture = integer(), activation = "sigmoid",
    decay = 0.001, loss = "logistic", restarts = 2,
    maxit = 200, seed = 17
  )
  expect_s3_class(affine, "odrl_nnet_backend_fit")
  expect_equal(affine$architecture, integer())
  expect_equal(affine$case_weight_sum, 1, tolerance = 1e-12)
  expect_true(is.finite(affine$objective))
  expect_true(all(is.finite(affine$all_objectives)))
  affine_score <- odrlITR:::.odrl_predict_nnet_backend(affine, x)
  expect_length(affine_score, nrow(x))
  expect_true(all(is.finite(affine_score)))
  expect_gt(mean(sign(affine_score) == score), 0.9)

  shallow <- odrlITR:::.odrl_fit_nnet_backend(
    x, score, architecture = 4L, activation = "sigmoid",
    decay = 0.001, loss = "logistic", restarts = 2,
    maxit = 200, seed = 18
  )
  expect_equal(shallow$architecture, 4L)
  expect_true(all(is.finite(
    odrlITR:::.odrl_predict_nnet_backend(shallow, x)
  )))
})

test_that("nnet entropy is the normalized weighted logistic margin risk", {
  skip_if_not_installed("nnet")
  set.seed(73)
  x <- matrix(rnorm(160), ncol = 2)
  score <- (0.5 + abs(rnorm(nrow(x)))) *
    sample(c(-1, 1), nrow(x), replace = TRUE)
  fit <- odrlITR:::.odrl_fit_nnet_backend(
    x, score, architecture = integer(), activation = "sigmoid",
    decay = 0, loss = "logistic", restarts = 1,
    maxit = 150, seed = 19
  )
  probability <- drop(stats::predict(fit$fit, x, type = "raw"))
  label <- as.numeric(score >= 0)
  case_weight <- abs(score) / (mean(abs(score)) * nrow(x))
  manual <- -sum(case_weight * (
    label * log(probability) + (1 - label) * log1p(-probability)
  ))
  expect_equal(fit$objective, manual, tolerance = 1e-10)
})

test_that("nnet backend restarts and score normalization are reproducible", {
  skip_if_not_installed("nnet")
  set.seed(72)
  x <- matrix(rnorm(120), ncol = 2)
  score <- ifelse(x[, 1] + x[, 2]^2 > 0.5, 2, -1)
  arguments <- list(
    x = x, score = score, architecture = 3L, activation = "sigmoid",
    decay = 0.01, loss = "logistic", restarts = 2,
    maxit = 150, seed = 91
  )
  first <- do.call(odrlITR:::.odrl_fit_nnet_backend, arguments)
  second <- do.call(odrlITR:::.odrl_fit_nnet_backend, arguments)
  scaled <- do.call(
    odrlITR:::.odrl_fit_nnet_backend,
    utils::modifyList(arguments, list(score = 100 * score))
  )

  expect_equal(first$fit$wts, second$fit$wts, tolerance = 0)
  expect_equal(first$all_objectives, second$all_objectives, tolerance = 0)
  expect_equal(first$fit$wts, scaled$fit$wts, tolerance = 1e-12)
  expect_equal(first$all_objectives, scaled$all_objectives,
               tolerance = 1e-12)
})

test_that("nnet backend factory follows the existing callback contract", {
  skip_if_not_installed("nnet")
  backend <- odrlITR:::.odrl_nnet_backend(list(skip = TRUE))
  expect_identical(backend$name, "nnet")
  expect_true(is.function(backend$fit))
  expect_true(is.function(backend$predict))

  x <- matrix(seq(-1, 1, length.out = 30), ncol = 1)
  score <- ifelse(x[, 1] >= 0, 1, -1)
  fit <- backend$fit(
    x = x, score = score, architecture = integer(),
    activation = "sigmoid", decay = 0.001, loss = "logistic",
    restarts = 1, maxit = 150, seed = 4, leaky_slope = 0.01
  )
  prediction <- backend$predict(fit, x)
  expect_length(prediction, nrow(x))
  expect_true(all(is.finite(prediction)))
})
