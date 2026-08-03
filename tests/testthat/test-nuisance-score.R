test_that("user nuisance produces the paper's score", {
  a <- c(-1, 1, -1, 1)
  y <- c(1, 2, 4, 8)
  m <- c(0.5, 1.5, 3, 7)
  pi <- c(0.25, 0.75, 0.4, 0.6)
  nuisance <- odrl_nuisance_user(
    m = m, pi = pi, out_of_fold = TRUE, source = "test truth"
  )
  observed <- odrl_score(a, y, nuisance)
  expected <- (a - (2 * pi - 1)) * (y - m)
  expect_equal(as.numeric(observed), expected)
  expect_equal(attr(observed, "working"), expected)
})

test_that("propensity clipping is explicit and audited", {
  nuisance <- odrl_nuisance_user(
    m = 1:3, pi = c(0, 0.5, 1), out_of_fold = TRUE,
    propensity_bounds = c(0.1, 0.9)
  )
  expect_equal(nuisance$pi, c(0.1, 0.5, 0.9))
  expect_equal(nuisance$diagnostics$clipped_fraction, 2 / 3)
})

test_that("scalar known propensity is recycled", {
  nuisance <- odrl_nuisance_user(
    m = 1:4, pi = 0.5, out_of_fold = TRUE, source = "randomized"
  )
  expect_equal(nuisance$pi, rep(0.5, 4))
})

test_that("factor positive level is explicit", {
  a <- factor(c("treated", "control", "treated", "control"))
  nuisance <- odrl_nuisance_user(
    m = rep(0, 4), pi = rep(0.5, 4), out_of_fold = TRUE
  )
  score <- odrl_score(a, rep(1, 4), nuisance, positive = "treated")
  expect_equal(as.numeric(score), c(1, -1, 1, -1))
})

test_that("nuisance metadata validation is fail closed", {
  expect_error(
    odrl_nuisance_user(m = 1:2, pi = c(0.4, 0.6), out_of_fold = NA),
    "TRUE or FALSE"
  )
  expect_error(
    odrl_nuisance_user(m = 1:2, pi = c(0.4, 0.6), source = character()),
    "nonempty"
  )
  expect_error(
    odrl_control(linear_require_gap = NA),
    "TRUE or FALSE"
  )
  expect_error(
    odrl_nuisance_user(m = factor(c("100", "200")), pi = 0.5),
    "genuinely numeric"
  )
  expect_error(
    odrl_nuisance_user(m = 1:2, pi = factor(c("0.2", "0.8"))),
    "genuinely numeric"
  )
  expect_error(
    odrl_nuisance_user(m = 1:2, pi = 0.5, fold_id = c(1, NA)),
    "cannot contain"
  )
})

test_that("a sparse score above the pointwise tolerance is not degenerate", {
  n <- 20
  a <- rep(c(-1, 1), n / 2)
  y <- rep(0, n)
  y[[1L]] <- -1e-9
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = 0.5, out_of_fold = TRUE
  )
  score <- odrl_score(a, y, nuisance, tolerance = 1e-10)
  expect_false(attr(score, "degenerate"))
  expect_equal(sum(attr(score, "working") != 0), 1)
})

test_that("score rejects a nuisance object with mismatched treatment coding", {
  a <- factor(rep(c("control", "treated"), 4))
  nuisance <- odrl_nuisance_user(
    m = rep(0, 8), pi = 0.5, out_of_fold = TRUE
  )
  nuisance$treatment_map <- odrl:::.odrl_encode_treatment(
    a, positive = "treated"
  )
  expect_error(
    odrl_score(a, rep(1, 8), nuisance, positive = "control"),
    "different treatment coding"
  )
})

test_that("score tolerance and numeric outcomes are validated", {
  nuisance <- odrl_nuisance_user(
    m = rep(0, 4), pi = 0.5, out_of_fold = TRUE
  )
  for (bad in list(NA_real_, -1, 0, Inf, c(1e-10, 1e-9))) {
    expect_error(odrl_score(c(-1, 1, -1, 1), 1:4, nuisance,
                            tolerance = bad), "tolerance")
  }
  expect_error(
    odrl_score(c(-1, 1, -1, 1), factor(1:4), nuisance),
    "numeric outcome"
  )
  expect_error(
    odrl(matrix(1:8, ncol = 2), c(-1, 1, -1, 1), factor(1:4),
         learner = "tree", nuisance = nuisance),
    "numeric outcome"
  )
})

test_that("propensity bounds cannot silently override a nuisance object", {
  nuisance <- odrl_nuisance_user(
    m = rep(0, 4), pi = c(0, 0.2, 0.8, 1), out_of_fold = TRUE
  )
  expect_error(
    odrl(matrix(1:8, ncol = 2), c(-1, 1, -1, 1), 1:4,
         learner = "tree", nuisance = nuisance,
         propensity_bounds = c(0.1, 0.9)),
    "cannot silently override"
  )
})

test_that("known propensity cannot silently override supplied nuisances", {
  nuisance <- odrl_nuisance_user(
    m = rep(0, 4), pi = 0.2, out_of_fold = TRUE
  )
  expect_error(
    odrl(matrix(1:8, ncol = 2), c(-1, 1, -1, 1), 1:4,
         learner = "tree", nuisance = nuisance, known_pi = 0.8),
    "only used with built-in nuisance"
  )
})
