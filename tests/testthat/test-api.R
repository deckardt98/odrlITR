test_that("zero scores use the deterministic positive tie rule", {
  x <- matrix(seq_len(20), ncol = 2)
  a <- rep(c(-1, 1), 5)
  y <- rep(3, 10)
  nuisance <- odrl_nuisance_user(
    m = y, pi = rep(0.5, 10), out_of_fold = TRUE
  )
  fit <- odrl(x, a, y, learner = "tree", nuisance = nuisance)
  expect_true(fit$score$degenerate)
  expect_equal(predict(fit, x), rep(1, 10))
})

test_that("unmarked user predictions trigger a cross-fitting warning", {
  dat <- odrl_simulate(60, seed = 3)
  nuisance <- odrl_nuisance_user(m = dat$m, pi = dat$pi)
  expect_warning(
    odrl(dat$x, dat$a, dat$y, learner = "relu", nuisance = nuisance,
         control = odrl_control(
           relu_hidden_units = 0, relu_decay = 0.01, relu_folds = 2,
           relu_restarts = 1, relu_refit_restarts = 1, relu_maxit = 150
         )),
    "not marked out-of-fold"
  )
})

test_that("formula interface excludes treatment from dot expansion", {
  dat <- odrl_simulate(60, seed = 8)
  frame <- data.frame(y = dat$y, a = dat$a, dat$x)
  nuisance <- odrl_nuisance_user(
    m = dat$m, pi = dat$pi, out_of_fold = TRUE
  )
  fit <- odrl_formula(
    y ~ ., treatment = "a", data = frame, learner = "relu",
    nuisance = nuisance,
    control = odrl_control(
      relu_hidden_units = 0, relu_decay = 0.01, relu_folds = 2,
      relu_restarts = 1, relu_refit_restarts = 1, relu_maxit = 150
    )
  )
  expect_equal(fit$p, ncol(dat$x))
})

test_that("formula predictions reproduce factors and interactions", {
  set.seed(22)
  n <- 80
  frame <- data.frame(
    y = rnorm(n),
    a = rep(c(-1, 1), n / 2),
    x = rnorm(n),
    g = factor(rep(c("a", "b"), each = n / 2))
  )
  nuisance <- odrl_nuisance_user(
    m = rep(0, n), pi = rep(0.5, n), out_of_fold = TRUE
  )
  fit <- odrl_formula(
    y ~ x * g, treatment = "a", data = frame, learner = "relu",
    nuisance = nuisance,
    control = odrl_control(
      relu_hidden_units = 0, relu_decay = 0.01, relu_folds = 2,
      relu_restarts = 1, relu_refit_restarts = 1, relu_maxit = 150
    )
  )
  newdata <- data.frame(x = c(-1, 1), g = factor(c("a", "b")))
  expect_length(predict(fit, newdata), 2)
  expect_equal(
    predict(fit, newdata),
    predict(fit, newdata[c("g", "x")])
  )
  expect_error(
    predict(fit, data.frame(x = 0, g = factor("unseen"))),
    "new level|new levels"
  )
  expect_equal(fit$p, 3)
})

test_that("named matrix predictions are invariant to column order", {
  x <- cbind(first = seq_len(10), second = seq_len(10)^2)
  a <- rep(c(-1, 1), 5)
  y <- rep(4, 10)
  nuisance <- odrl_nuisance_user(
    m = y, pi = 0.5, out_of_fold = TRUE
  )
  fit <- odrl(x, a, y, learner = "tree", nuisance = nuisance)
  expect_equal(predict(fit, x), predict(fit, x[, 2:1, drop = FALSE]))
  expect_error(
    predict(fit, cbind(first = x[, 1], wrong = x[, 2])),
    "do not match training columns"
  )
})

test_that("character treatments are returned as character actions", {
  x <- matrix(seq_len(20), ncol = 2)
  a <- rep(c("control", "treated"), 5)
  y <- rep(2, 10)
  nuisance <- odrl_nuisance_user(
    m = y, pi = 0.5, out_of_fold = TRUE
  )
  fit <- odrl(x, a, y, learner = "tree", nuisance = nuisance,
              positive = "treated")
  action <- predict(fit, x)
  expect_type(action, "character")
  expect_equal(action, rep("treated", 10))
})

test_that("formula treatment column name is scalar and nonempty", {
  frame <- data.frame(y = 1:4, a = c(-1, 1, -1, 1), x = 1:4)
  expect_error(odrl_formula(y ~ x, character(), frame), "one nonempty")
  expect_error(odrl_formula(y ~ x, c("a", "x"), frame), "one nonempty")
})
