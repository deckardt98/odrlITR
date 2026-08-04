test_that("series specifications expose conservative aliases and grids", {
  expect_error(odrl_series_kernel(character()), "basis")
  spline <- odrlITR:::.odrl_resolve_series_kernel("spline")
  wavelet <- odrlITR:::.odrl_resolve_series_kernel("wavelet")
  local <- odrlITR:::.odrl_resolve_series_kernel("local-poly")

  expect_s3_class(spline, "odrl_series_kernel")
  expect_identical(spline$basis, "bspline")
  expect_identical(wavelet$basis, "haar")
  expect_identical(local$basis, "local_polynomial")
  expect_equal(nrow(odrlITR:::.odrl_series_grid(spline)), 3L)
  expect_equal(nrow(odrlITR:::.odrl_series_grid(local)), 4L)

  expect_error(
    odrlITR:::odrl_series_kernel(
      basis = "fourier", combine = "total_degree"
    ),
    "requires a Legendre"
  )
  expect_error(
    odrlITR:::odrl_series_kernel(boundary_quantiles = c(0.9, 0.1)),
    "increasing"
  )
})

test_that("univariate series agree with hand-calculable values", {
  u <- c(0, 0.5, 1)
  legendre <- odrlITR:::.odrl_legendre_matrix(u, 2L)
  expect_equal(legendre[, 1L], sqrt(3) * c(-1, 0, 1))
  expect_equal(legendre[, 2L], sqrt(5) * c(1, -0.5, 1))

  fourier <- odrlITR:::.odrl_fourier_matrix(c(0, 0.25), 1L)
  expect_equal(fourier[, 1L], sqrt(2) * c(1, 0), tolerance = 1e-12)
  expect_equal(fourier[, 2L], sqrt(2) * c(0, 1), tolerance = 1e-12)

  haar <- odrlITR:::.odrl_haar_matrix(c(0.25, 0.75), 1L)
  expect_equal(drop(haar), c(1, -1))

  local <- odrlITR:::.odrl_local_polynomial_matrix(
    c(0.25, 0.75), partitions = 2L, degree = 1L
  )
  expect_equal(dim(local), c(2L, 3L))
  expect_equal(local[, 1L], c(1, 0))
})

test_that("fitted maps persist centering, RMS scaling, and a PSD feature kernel", {
  set.seed(801)
  x <- matrix(runif(180), ncol = 3)
  spec <- odrlITR:::odrl_series_kernel(
    basis = "legendre", legendre_degree = 2L,
    domain = "unit", bounds = c(0, 1), max_features = 100L
  )
  mapped <- odrlITR:::.odrl_fit_series_map(x, spec, list(degree = 2L))
  features <- mapped$features

  expect_equal(mapped$map$raw_feature_count, 6L)
  expect_equal(mapped$map$feature_count, 6L)
  expect_equal(unname(colMeans(features)), rep(0, 6), tolerance = 1e-12)
  expect_equal(
    unname(colMeans(features^2)), rep(1 / 6, 6), tolerance = 1e-10
  )
  kernel <- tcrossprod(features)
  expect_equal(kernel, t(kernel), tolerance = 1e-12)
  expect_gte(min(eigen(kernel, symmetric = TRUE, only.values = TRUE)$values),
             -1e-9)
  expect_equal(
    odrlITR:::.odrl_apply_series_map(mapped$map, x), features,
    tolerance = 1e-12
  )
})

test_that("multivariate construction counts are exact", {
  set.seed(802)
  x <- matrix(runif(240), ncol = 3)
  fit_count <- function(combine, order = 1L) {
    spec <- odrlITR:::odrl_series_kernel(
      "legendre", legendre_degree = 2L, combine = combine,
      interaction_order = order, bounds = c(0, 1), domain = "unit",
      max_features = 100L
    )
    odrlITR:::.odrl_fit_series_map(
      x, spec, list(degree = 2L)
    )$map$raw_feature_count
  }

  expect_equal(fit_count("additive"), 6L)
  expect_equal(fit_count("anova", 2L), 18L)
  expect_equal(fit_count("total_degree"), 9L)
  expect_equal(fit_count("tensor"), 26L)
})

test_that("feature explosion is rejected before allocation", {
  x <- matrix(runif(200), ncol = 10)
  spec <- odrlITR:::odrl_series_kernel(
    "legendre", legendre_degree = 3L, combine = "tensor",
    bounds = c(0, 1), domain = "unit", max_features = 1000L
  )
  expect_error(
    odrlITR:::.odrl_fit_series_map(x, spec, list(degree = 3L)),
    "exceeding `max_features"
  )

  elements <- odrlITR:::odrl_series_kernel(
    "fourier", fourier_harmonics = 2L, bounds = c(0, 1),
    domain = "unit", max_features = 100L, max_feature_elements = 10L
  )
  expect_error(
    odrlITR:::.odrl_fit_series_map(
      matrix(runif(40), ncol = 2), elements, list(harmonics = 2L)
    ),
    "max_feature_elements"
  )
})

test_that("training-fold bounds and spline knots are retained for prediction", {
  x <- cbind(x1 = seq(0, 1, length.out = 40), x2 = runif(40))
  spec <- odrlITR:::odrl_series_kernel(
    "bspline", spline_df = 5L, spline_degree = 3L,
    bounds = c(0, 1), domain = "unit"
  )
  mapped <- odrlITR:::.odrl_fit_series_map(x, spec, list(df = 5L))
  knots <- mapped$map$univariate_state[[1L]]$knots
  prediction <- odrlITR:::.odrl_apply_series_map(mapped$map, x[1:5, ])

  expect_equal(
    mapped$map$univariate_state[[1L]]$knots, knots, tolerance = 0
  )
  expect_true(all(is.finite(prediction)))
  expect_equal(ncol(prediction), mapped$map$feature_count)

  clamped <- odrlITR:::.odrl_apply_series_map(
    mapped$map, rbind(c(-10, 0.5), c(10, 0.5)), diagnostics = TRUE
  )
  expect_gt(sum(clamped$clamped_below), 0)
  expect_gt(sum(clamped$clamped_above), 0)

  error_spec <- odrlITR:::odrl_series_kernel(
    "bspline", spline_df = 5L, bounds = c(0, 1), domain = "unit",
    extrapolation = "error"
  )
  error_map <- odrlITR:::.odrl_fit_series_map(
    x, error_spec, list(df = 5L)
  )$map
  expect_error(
    odrlITR:::.odrl_apply_series_map(error_map, rbind(c(-1, 0.5))),
    "outside"
  )
})

test_that("all built-in series fit finite primal surrogate rules", {
  set.seed(803)
  x <- matrix(runif(180), ncol = 3)
  score <- x[, 1L] - x[, 2L] + 0.2 * sin(2 * pi * x[, 3L])
  specifications <- list(
    list(
      spec = odrlITR:::odrl_series_kernel(
        "legendre", legendre_degree = 2L, bounds = c(0, 1)
      ),
      candidate = list(degree = 2L)
    ),
    list(
      spec = odrlITR:::odrl_series_kernel(
        "fourier", fourier_harmonics = 1L, bounds = c(0, 1)
      ),
      candidate = list(harmonics = 1L)
    ),
    list(
      spec = odrlITR:::odrl_series_kernel(
        "bspline", spline_df = 4L, bounds = c(0, 1)
      ),
      candidate = list(df = 4L)
    ),
    list(
      spec = odrlITR:::odrl_series_kernel(
        "haar", wavelet_level = 2L, bounds = c(0, 1)
      ),
      candidate = list(level = 2L)
    ),
    list(
      spec = odrlITR:::odrl_series_kernel(
        "local_polynomial", local_partitions = 2L, local_degree = 1L,
        bounds = c(0, 1)
      ),
      candidate = list(partitions = 2L, local_degree = 1L)
    )
  )

  for (item in specifications) {
    fit <- odrlITR:::.odrl_fit_series_surrogate(
      x, score, item$spec, item$candidate, lambda = 0.1,
      maxit = 500L, seed = 11L, loss = "logistic"
    )
    prediction <- odrlITR:::.odrl_predict_series_unclipped(fit, x)
    expect_true(is.finite(fit$objective), info = item$spec$basis)
    expect_true(all(is.finite(prediction)), info = item$spec$basis)
    expect_equal(length(prediction), nrow(x), info = item$spec$basis)
  }
})

test_that("series maps survive serialization without prediction drift", {
  set.seed(804)
  x <- matrix(runif(120), ncol = 2)
  score <- x[, 1L] * x[, 2L] - 0.2
  spec <- odrlITR:::odrl_series_kernel(
    "legendre", legendre_degree = 2L, combine = "anova",
    interaction_order = 2L, bounds = c(0, 1)
  )
  fit <- odrlITR:::.odrl_fit_series_surrogate(
    x, score, spec, list(degree = 2L), lambda = 0.1,
    maxit = 500L, seed = 12L, loss = "exponential"
  )
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(
    odrlITR:::.odrl_predict_series_unclipped(restored, x),
    odrlITR:::.odrl_predict_series_unclipped(fit, x),
    tolerance = 0
  )
})
