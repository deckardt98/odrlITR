test_that("Super Learner nuisances are outer-fold predictions", {
  skip_if_not_installed("SuperLearner")
  dat <- odrl_simulate(80, seed = 15)
  nuisance <- odrl_nuisance_sl(
    dat$x, dat$a, dat$y,
    folds = 2, inner_folds = 2,
    sl.library = c("SL.mean", "SL.glm"), seed = 2
  )
  expect_s3_class(nuisance, "odrl_nuisance")
  expect_true(nuisance$out_of_fold)
  expect_setequal(unique(nuisance$fold_id), 1:2)
  expect_length(nuisance$m, 80)
  expect_true(all(nuisance$pi >= 0 & nuisance$pi <= 1))
})

test_that("known propensity can be combined with cross-fitted m", {
  skip_if_not_installed("SuperLearner")
  dat <- odrl_simulate(80, seed = 44)
  nuisance <- odrl_nuisance_sl(
    dat$x, dat$a, dat$y, folds = 2, inner_folds = 2,
    sl.library.m = c("SL.mean", "SL.glm"), known_pi = 0.5, seed = 3
  )
  expect_true(nuisance$known_propensity)
  expect_equal(nuisance$pi, rep(0.5, 80))
  expect_true(nuisance$out_of_fold)
})

test_that("inner propensity folds respect arm counts", {
  skip_if_not_installed("SuperLearner")
  set.seed(30)
  n <- 100
  x <- matrix(rnorm(n * 2), ncol = 2)
  a <- c(rep(1, 10), rep(-1, 90))
  y <- rnorm(n)
  nuisance <- odrl_nuisance_sl(
    x, a, y, folds = 2, inner_folds = 10,
    sl.library = c("SL.mean", "SL.glm"), seed = 5
  )
  expect_s3_class(nuisance, "odrl_nuisance")
  expect_true(all(vapply(nuisance$fits, function(x) x$holdout_n > 0,
                         logical(1))))
})

test_that("known propensity permits imbalanced arms in pooled m fitting", {
  skip_if_not_installed("SuperLearner")
  set.seed(61)
  n <- 40
  x <- matrix(rnorm(n * 2), ncol = 2)
  a <- c(1, 1, rep(-1, n - 2))
  y <- rnorm(n)
  nuisance <- odrl_nuisance_sl(
    x, a, y, folds = 2, inner_folds = 3,
    sl.library.m = "SL.mean", known_pi = 0.05, seed = 8
  )
  expect_true(nuisance$known_propensity)
  expect_equal(nuisance$pi, rep(0.05, n))
  expect_true(all(vapply(nuisance$fits, function(z) {
    z$outcome_inner_folds == 3L
  }, logical(1))))
})

test_that("captured Super Learner warnings remain auditable", {
  captured <- odrlITR:::.odrl_capture_warnings({
    warning("deliberate test warning")
    7
  })
  expect_equal(captured$value, 7)
  expect_match(captured$warnings, "deliberate test warning")
})
