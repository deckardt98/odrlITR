.small_control <- function() {
  odrl_control(
    svm_rbf_multiplier = 1,
    svm_penalty = 0.1,
    svm_folds = 2,
    svm_maxit = 150,
    relu_hidden_units = c(0, 3),
    relu_decay = 0.01,
    relu_folds = 2,
    relu_restarts = 1,
    relu_refit_restarts = 1,
    relu_maxit = 150,
    seed = 9
  )
}

test_that("bounded Gaussian hinge scores use hard-tanh clipping", {
  dat <- odrl_simulate(60, seed = 2)
  nuisance <- odrl_nuisance_user(
    m = dat$m, pi = dat$pi, out_of_fold = TRUE
  )
  fit <- odrl(
    dat$x, dat$a, dat$y, learner = "svm", loss = "hinge",
    nuisance = nuisance, control = .small_control()
  )
  newx <- matrix(rnorm(200), ncol = ncol(dat$x))
  colnames(newx) <- colnames(dat$x)
  score <- predict(fit, newx, type = "score")
  expect_lte(max(abs(score)), 1 + 1e-8)
  expect_true(fit$policy$diagnostics$globally_bounded)
  z <- odrlITR:::.odrl_standardize_apply(newx, fit$policy$transform)
  raw <- odrlITR:::.odrl_predict_kernel_unclipped(
    fit$policy$fit, z, fit$policy$kernel
  )
  expect_equal(score, odrlITR:::.odrl_hardtanh(raw), tolerance = 1e-8)
  expect_identical(fit$policy$diagnostics$hinge_mode, "bounded")
  expect_identical(fit$policy$diagnostics$clipping, "hard_tanh")
  expect_true(is.finite(fit$policy$fit$rkhs_norm))
  expect_true("lambda" %in% names(fit$policy$selected))
})

test_that("kernel logistic returns finite scores", {
  dat <- odrl_simulate(50, seed = 4)
  nuisance <- odrl_nuisance_user(
    m = dat$m, pi = dat$pi, out_of_fold = TRUE
  )
  fit <- odrl(
    dat$x, dat$a, dat$y, learner = "svm", loss = "logistic",
    nuisance = nuisance, control = .small_control()
  )
  expect_true(all(is.finite(predict(fit, dat$x, type = "score"))))
})

test_that("ReLU supports bounded hinge and logistic losses", {
  dat <- odrl_simulate(60, seed = 6)
  nuisance <- odrl_nuisance_user(
    m = dat$m, pi = dat$pi, out_of_fold = TRUE
  )
  hinge <- odrl(
    dat$x, dat$a, dat$y, learner = "relu", loss = "hinge",
    nuisance = nuisance, control = .small_control()
  )
  logistic <- odrl(
    dat$x, dat$a, dat$y, learner = "relu", loss = "logistic",
    nuisance = nuisance, control = .small_control()
  )
  expect_lte(max(abs(predict(hinge, dat$x, type = "score"))), 1)
  expect_true(all(is.finite(predict(logistic, dat$x, type = "score"))))
})

test_that("hinge ReLU initializes inside its active region", {
  x <- matrix(c(rep(-2, 10), rep(2, 90)), ncol = 1)
  score <- c(rep(-1, 10), rep(1, 90))
  fit <- odrlITR:::.odrl_fit_relu_candidate(
    x, score, hidden = 0, decay = 0.001, loss = "hinge",
    restarts = 1, maxit = 200, seed = 19
  )
  action <- ifelse(odrlITR:::.odrl_predict_relu_raw(fit, x) >= 0, 1, -1)
  expect_equal(fit$convergence, 0L)
  expect_lte(mean(action != sign(score)), 0.01)
})
