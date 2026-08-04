test_that("built-in folds are treatment-stratified and size-balanced", {
  a <- c(rep(-1, 53), rep(1, 48))
  fold_id <- odrlITR:::.odrl_treatment_stratified_folds(a, 5, seed = 18)

  total <- tabulate(fold_id, nbins = 5)
  negative <- tabulate(fold_id[a == -1], nbins = 5)
  positive <- tabulate(fold_id[a == 1], nbins = 5)

  expect_lte(max(total) - min(total), 1)
  expect_lte(max(negative) - min(negative), 1)
  expect_lte(max(positive) - min(positive), 1)
  expect_true(all(negative > 0))
  expect_true(all(positive > 0))
  expect_identical(
    fold_id,
    odrlITR:::.odrl_treatment_stratified_folds(a, 5, seed = 18)
  )
})

test_that("parametric nuisances cross-fit logistic and linear models", {
  set.seed(42)
  n <- 180
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  pi <- plogis(-0.2 + 0.7 * x[, 1] - 0.4 * x[, 2])
  a <- rbinom(n, 1, pi)
  y <- 1.5 + 1.2 * x[, 1] - 0.8 * x[, 2] + rnorm(n, sd = 0.1)

  nuisance <- odrlITR:::odrl_nuisance_parametric(
    x, a, y, folds = 3, seed = 6
  )

  expect_s3_class(nuisance, "odrl_nuisance")
  expect_true(nuisance$out_of_fold)
  expect_identical(nuisance$engine, "parametric")
  expect_identical(nuisance$folds, 3L)
  expect_identical(nuisance$fold_source, "treatment-stratified")
  expect_length(nuisance$fits, 3)
  expect_true(all(is.finite(nuisance$m)))
  expect_true(all(nuisance$pi >= 0 & nuisance$pi <= 1))
  expect_gt(stats::cor(nuisance$m, y), 0.98)
  expect_true(all(vapply(nuisance$fits, function(fit) {
    identical(fit$propensity_engine, "binomial logistic regression") &&
      identical(fit$outcome_engine, "Gaussian linear regression")
  }, logical(1))))
})

test_that("user-specified folds are respected and audited", {
  set.seed(91)
  n <- 120
  x <- matrix(rnorm(n * 2), ncol = 2)
  a <- rep(c(-1, 1), n / 2)
  y <- 1 + x[, 1] + rnorm(n)
  requested <- factor(rep(c("alpha", "beta", "gamma"), length.out = n),
                      levels = c("gamma", "beta", "alpha"))

  nuisance <- odrlITR:::odrl_nuisance_parametric(
    x, a, y, folds = 17, fold_id = requested, seed = 3
  )

  expect_identical(nuisance$folds, 3L)
  expect_identical(nuisance$fold_source, "user supplied")
  expect_identical(nuisance$fold_labels, c("alpha", "beta", "gamma"))
  expect_identical(nuisance$fold_id, match(as.character(requested),
                                           c("alpha", "beta", "gamma")))
  expect_identical(
    vapply(nuisance$fits, `[[`, character(1), "fold_label"),
    c("alpha", "beta", "gamma")
  )
})

test_that("user-specified fold identifiers fail closed", {
  x <- matrix(seq_len(40), ncol = 2)
  a <- rep(c(-1, 1), 10)
  y <- seq_len(20)

  expect_error(
    odrlITR:::odrl_nuisance_parametric(
      x, a, y, fold_id = rep(1, 19)
    ),
    "must match the number of rows"
  )
  expect_error(
    odrlITR:::odrl_nuisance_parametric(
      x, a, y, fold_id = rep(1, 20)
    ),
    "at least two"
  )
  expect_error(
    odrlITR:::odrl_nuisance_parametric(
      x, a, y, fold_id = c(rep("a", 19), NA_character_)
    ),
    "cannot contain missing"
  )

  # The first training set contains only controls and the second only treated.
  separated_a <- c(rep(-1, 10), rep(1, 10))
  separated_fold <- c(rep("left", 10), rep("right", 10))
  expect_error(
    odrlITR:::odrl_nuisance_parametric(
      x, separated_a, y, fold_id = separated_fold
    ),
    "training set must contain at least two"
  )
})

test_that("known propensity works with parametric cross-fitted outcome model", {
  set.seed(71)
  n <- 60
  x <- matrix(rnorm(n * 2), ncol = 2)
  a <- c(1, 1, rep(-1, n - 2))
  y <- 2 + x[, 1] + rnorm(n)
  fold_id <- rep(c("one", "two"), length.out = n)

  nuisance <- odrlITR:::odrl_nuisance_parametric(
    x, a, y, fold_id = fold_id, known_pi = 0.05, seed = 2
  )

  expect_true(nuisance$known_propensity)
  expect_equal(nuisance$pi, rep(0.05, n))
  expect_true(nuisance$out_of_fold)

  # With a known propensity, automatic folds need not place the single
  # treated observation in every outer training sample: only pooled m is fit.
  highly_imbalanced <- odrlITR:::odrl_nuisance_parametric(
    x, c(1, rep(-1, n - 1)), y, folds = 5, known_pi = 0.05, seed = 2
  )
  expect_identical(highly_imbalanced$folds, 5L)
  expect_lte(max(tabulate(highly_imbalanced$fold_id)) -
               min(tabulate(highly_imbalanced$fold_id)), 1)
})

test_that("built-in nuisance fitting rejects undersized outer training sets", {
  x <- matrix(seq_len(10), ncol = 2)
  a <- c(1, rep(-1, 4))
  y <- seq_len(5)
  pathological <- c("singleton", rep("remainder", 4))

  expect_error(
    odrlITR:::odrl_nuisance_parametric(
      x, a, y, fold_id = pathological, known_pi = 0.2
    ),
    "at least two observations"
  )
  skip_if_not_installed("SuperLearner")
  expect_error(
    odrl_nuisance_sl(
      x, a, y, fold_id = pathological, known_pi = 0.2,
      inner_folds = 2, sl.library = "SL.mean"
    ),
    "at least two observations"
  )
})

test_that("Super Learner accepts user-specified outer folds", {
  skip_if_not_installed("SuperLearner")
  dat <- odrl_simulate(80, seed = 37)
  fold_id <- rep(c("A", "B"), length.out = 80)
  nuisance <- odrl_nuisance_sl(
    dat$x, dat$a, dat$y,
    folds = 9, fold_id = fold_id, inner_folds = 2,
    sl.library = c("SL.mean", "SL.glm"), seed = 2
  )

  expect_identical(nuisance$folds, 2L)
  expect_identical(nuisance$fold_source, "user supplied")
  expect_identical(nuisance$fold_labels, c("A", "B"))
  expect_identical(nuisance$fold_id, rep(1:2, length.out = 80))
})

test_that("user nuisance fold metadata is strictly validated", {
  expect_error(
    odrl_nuisance_user(m = 1:4, pi = 0.5, fold_id = rep("only", 4)),
    "at least two"
  )
  expect_error(
    odrl_nuisance_user(m = 1:4, pi = 0.5,
                       fold_id = c("a", "b", "", "b")),
    "empty labels"
  )
})
