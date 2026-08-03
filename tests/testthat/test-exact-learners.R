test_that("direct affine MILP passes sign and objective audit", {
  skip_if_not_installed("highs")
  x <- matrix(c(-2, -1, 1, 2), ncol = 1,
              dimnames = list(NULL, "x"))
  score <- c(-2, -1, 1, 2)
  fit <- odrl:::.odrl_fit_linear(
    x, score,
    odrl_control(linear_time_limit = 10, linear_relative_gap = 0)
  )
  expect_equal(fit$training_action, c(-1, -1, 1, 1))
  expect_equal(fit$diagnostics$action_mismatches, 0)
  expect_lte(fit$diagnostics$max_constraint_violation, 1e-7)
  expect_lte(fit$diagnostics$objective_difference, 1e-7)
  expect_identical(
    fit$diagnostics$proved_optimal,
    isTRUE(fit$diagnostics$solver_reported_optimal) &&
      is.finite(fit$diagnostics$mip_gap) &&
      fit$diagnostics$mip_gap <= 1e-8
  )
  expect_lte(fit$diagnostics$max_integrality_error,
             fit$diagnostics$integrality_tolerance)
})

test_that("policy tree returns binary actions", {
  skip_if_not_installed("policytree")
  dat <- odrl_simulate(100, boundary = "tree", seed = 12)
  nuisance <- odrl_nuisance_user(
    m = dat$m, pi = dat$pi, out_of_fold = TRUE
  )
  fit <- odrl(
    dat$x, dat$a, dat$y, learner = "tree", nuisance = nuisance,
    control = odrl_control(tree_min_node_size = 2)
  )
  expect_true(all(predict(fit, dat$x) %in% c(-1, 1)))
  expect_equal(length(fitted(fit)), 100)
})

test_that("policy tree resolves tied leaf rewards to positive treatment", {
  skip_if_not_installed("policytree")
  x <- matrix(1:4, ncol = 1)
  score <- c(-2, -1, 1, 2)
  fit <- odrl:::.odrl_fit_tree(
    x, score,
    odrl_control(tree_depth = 0, tree_min_node_size = 1)
  )
  expect_equal(fit$training_action, rep(1, 4))
  expect_true(fit$diagnostics$global_candidate_search)
})

test_that("affine coefficients default to the original covariate scale", {
  skip_if_not_installed("highs")
  x <- cbind(x1 = seq(-20, 20, length.out = 20),
             x2 = seq(100, 300, length.out = 20)^2)
  a <- rep(c(-1, 1), 10)
  y <- a * x[, "x1"]
  nuisance <- odrl_nuisance_user(
    m = rep(0, 20), pi = 0.5, out_of_fold = TRUE
  )
  fit <- odrl(
    x, a, y, learner = "linear", nuisance = nuisance,
    control = odrl_control(linear_time_limit = 10)
  )
  beta <- coef(fit)
  expect_named(beta, c("(Intercept)", "x1", "x2"))
  reconstructed <- drop(cbind(1, x) %*% beta)
  expect_equal(reconstructed, predict(fit, x, type = "score"),
               tolerance = 1e-7)
  expect_equal(
    reconstructed,
    predict(fit, x[, c("x2", "x1")], type = "score"),
    tolerance = 1e-7
  )
  expect_equal(coef(fit, standardized = TRUE),
               fit$policy$coefficients)
  score <- predict(fit, x, type = "score")
  expect_equal(predict(fit, x), ifelse(score >= 0, 1, -1))
})

test_that("affine control rejects an infeasible margin convention", {
  expect_error(
    odrl_control(linear_coefficient_bound = 1, linear_margin = 2),
    "cannot exceed"
  )
})
