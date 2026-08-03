.odrl_squared_distance <- function(x, y = x) {
  x2 <- rowSums(x * x)
  y2 <- rowSums(y * y)
  pmax(outer(x2, y2, "+") - 2 * tcrossprod(x, y), 0)
}

.odrl_median_squared_distance <- function(x, seed, max_rows = 500L) {
  if (nrow(x) > max_rows) {
    set.seed(seed)
    x <- x[sample.int(nrow(x), max_rows), , drop = FALSE]
  }
  d <- .odrl_squared_distance(x)
  d <- d[upper.tri(d) & is.finite(d) & d > 1e-12]
  if (!length(d)) return(1)
  stats::median(d)
}

.odrl_kernel_matrix <- function(x, y = x, kernel, bandwidth2 = NULL) {
  if (kernel == "linear") {
    return(tcrossprod(x, y) / max(1, ncol(x)))
  }
  exp(-.odrl_squared_distance(x, y) / (2 * bandwidth2))
}

.odrl_fit_kernel_hinge <- function(x, score, kernel, multiplier, radius,
                                   seed) {
  if (kernel != "rbf") {
    .odrl_abort(
      "Bounded-hinge SVM currently requires `svm_kernel = \"rbf\"`. ",
      "For a Gaussian kernel with K(x,x)=1 and radius <= 1, the fitted ",
      "score has the paper's global [-1,1] certificate."
    )
  }
  bandwidth2 <- .odrl_median_squared_distance(x, seed) * multiplier
  k <- .odrl_kernel_matrix(x, kernel = kernel, bandwidth2 = bandwidth2)
  direction <- score / length(score)
  norm <- sqrt(max(drop(crossprod(direction, k %*% direction)), 0))
  alpha <- if (norm <= 1e-14) rep(0, length(score)) else {
    radius * direction / norm
  }
  fitted <- drop(k %*% alpha)
  bounded <- .odrl_hardtanh(fitted)
  list(
    alpha = alpha,
    training_x = x,
    bandwidth2 = bandwidth2,
    radius = radius,
    fitted = bounded,
    fitted_unclipped = fitted,
    convergence = 0L,
    objective = -mean(score * bounded),
    rkhs_norm = if (norm <= 1e-14) 0 else radius,
    globally_bounded = TRUE
  )
}

.odrl_fit_kernel_logistic <- function(x, score, kernel, multiplier, lambda,
                                      maxit, seed) {
  bandwidth2 <- if (kernel == "rbf") {
    .odrl_median_squared_distance(x, seed) * multiplier
  } else {
    NA_real_
  }
  k <- .odrl_kernel_matrix(x, kernel = kernel, bandwidth2 = bandwidth2)
  label <- ifelse(score >= 0, 1, -1)
  weight <- abs(score) / mean(abs(score))
  n <- length(score)
  objective <- function(theta) {
    intercept <- theta[[1L]]
    alpha <- theta[-1L]
    raw <- intercept + drop(k %*% alpha)
    mean(weight * .odrl_log1pexp(-label * raw)) +
      0.5 * lambda * (intercept^2 +
        drop(crossprod(alpha, k %*% alpha)))
  }
  gradient <- function(theta) {
    intercept <- theta[[1L]]
    alpha <- theta[-1L]
    raw <- intercept + drop(k %*% alpha)
    derivative <- -weight * label * stats::plogis(-label * raw) / n
    c(sum(derivative) + lambda * intercept,
      drop(k %*% derivative) + lambda * drop(k %*% alpha))
  }
  set.seed(seed)
  initial_intercept <- stats::qlogis(pmin(pmax(mean(label == 1), 0.01), 0.99))
  fit <- stats::optim(
    par = c(initial_intercept, rep(0, n)), fn = objective, gr = gradient,
    method = "L-BFGS-B", control = list(maxit = maxit, factr = 1e8)
  )
  attempts <- 1L
  if (fit$convergence != 0L) {
    retry <- stats::optim(
      par = fit$par, fn = objective, gr = gradient,
      method = "L-BFGS-B",
      control = list(maxit = max(3L * maxit, maxit + 100L), factr = 1e8)
    )
    attempts <- 2L
    if (retry$convergence == 0L || retry$value < fit$value) fit <- retry
  }
  list(
    intercept = fit$par[[1L]],
    alpha = fit$par[-1L],
    training_x = x,
    bandwidth2 = bandwidth2,
    lambda = lambda,
    fitted = fit$par[[1L]] + drop(k %*% fit$par[-1L]),
    convergence = fit$convergence,
    message = fit$message,
    attempts = attempts,
    objective = fit$value,
    globally_bounded = FALSE
  )
}

.odrl_predict_kernel_unclipped <- function(fit, newx, kernel) {
  k <- .odrl_kernel_matrix(
    newx, fit$training_x, kernel = kernel, bandwidth2 = fit$bandwidth2
  )
  (fit$intercept %||% 0) + drop(k %*% fit$alpha)
}

.odrl_predict_kernel_raw <- function(fit, newx, kernel, loss) {
  raw <- .odrl_predict_kernel_unclipped(fit, newx, kernel)
  if (loss == "hinge") .odrl_hardtanh(raw) else raw
}

.odrl_kernel_grid <- function(control, loss) {
  multiplier <- if (control$svm_kernel == "rbf") {
    control$svm_rbf_multiplier
  } else {
    1
  }
  if (loss == "hinge") {
    expand.grid(
      multiplier = multiplier,
      radius = control$svm_radius,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    expand.grid(
      multiplier = multiplier,
      lambda = rev(control$svm_penalty),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

.odrl_fit_svm_candidate <- function(x, score, control, loss, candidate, seed) {
  if (loss == "hinge") {
    .odrl_fit_kernel_hinge(
      x, score, control$svm_kernel, candidate$multiplier,
      candidate$radius, seed
    )
  } else {
    .odrl_fit_kernel_logistic(
      x, score, control$svm_kernel, candidate$multiplier,
      candidate$lambda, control$svm_maxit, seed
    )
  }
}

#' Fit a tuned kernel ODRL rule
#' @noRd
.odrl_fit_svm <- function(x, score, control, loss) {
  if (loss == "hinge" && control$svm_kernel != "rbf") {
    .odrl_abort("Bounded-hinge SVM requires `svm_kernel = \"rbf\"`.")
  }
  folds <- min(control$svm_folds, nrow(x))
  fold_id <- .odrl_score_folds(score, folds, control$seed + 300L)
  grid <- .odrl_kernel_grid(control, loss)
  criterion <- matrix(NA_real_, nrow(grid), folds)
  convergence <- matrix(NA_integer_, nrow(grid), folds)
  started <- proc.time()[["elapsed"]]
  for (j in seq_len(nrow(grid))) {
    candidate <- as.list(grid[j, , drop = FALSE])
    for (fold in seq_len(folds)) {
      train <- fold_id != fold
      holdout <- !train
      transform <- .odrl_standardize_fit(x[train, , drop = FALSE])
      train_x <- .odrl_standardize_apply(x[train, , drop = FALSE], transform)
      holdout_x <- .odrl_standardize_apply(x[holdout, , drop = FALSE], transform)
      fit <- .odrl_fit_svm_candidate(
        train_x, score[train], control, loss, candidate,
        control$seed + 10000L * j + fold
      )
      raw <- .odrl_predict_kernel_raw(
        fit, holdout_x, control$svm_kernel, loss
      )
      action <- ifelse(raw >= 0, 1, -1)
      criterion[j, fold] <- .odrl_empirical_criterion(score[holdout], action)
      convergence[j, fold] <- fit$convergence
    }
  }
  grid$mean_criterion <- rowMeans(criterion)
  grid$se_criterion <- apply(criterion, 1L, stats::sd) / sqrt(folds)
  grid$convergence_failures <- rowSums(convergence != 0L, na.rm = TRUE)
  eligible <- which(
    grid$convergence_failures == 0L & is.finite(grid$mean_criterion)
  )
  if (!length(eligible)) {
    .odrl_abort(
      "Every kernel candidate had at least one nonconverged tuning fit. ",
      "Increase `svm_maxit` or simplify the tuning grid."
    )
  }
  best <- eligible[[which.max(grid$mean_criterion[eligible])]]
  transform <- .odrl_standardize_fit(x)
  standardized_x <- .odrl_standardize_apply(x, transform)
  final <- .odrl_fit_svm_candidate(
    standardized_x, score, control, loss,
    as.list(grid[best, , drop = FALSE]), control$seed + 900000L
  )
  if (final$convergence != 0L) {
    .odrl_abort(
      "The selected kernel candidate did not converge on the full sample. ",
      "Increase `svm_maxit`."
    )
  }
  raw <- .odrl_predict_kernel_raw(
    final, standardized_x, control$svm_kernel, loss
  )
  action <- ifelse(raw >= 0, 1, -1)
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    engine = "kernel",
    kernel = control$svm_kernel,
    loss = loss,
    fit = final,
    transform = transform,
    tuning = grid,
    selected = grid[best, , drop = FALSE],
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      folds = folds,
      elapsed = elapsed,
      globally_bounded = isTRUE(final$globally_bounded),
      rkhs_norm = final$rkhs_norm %||% NA_real_,
      radius = final$radius %||% NA_real_,
      max_training_abs_unclipped = if (loss == "hinge") {
        max(abs(final$fitted_unclipped))
      } else {
        NA_real_
      },
      optimization_score_scale = "mean-absolute-score normalized",
      convergence = final$convergence
    )
  ), class = "odrl_policy_svm")
}

.odrl_predict_svm <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
  z <- .odrl_standardize_apply(newx, object$transform)
  raw <- .odrl_predict_kernel_raw(object$fit, z, object$kernel, object$loss)
  if (type == "score") raw else ifelse(raw >= 0, 1, -1)
}
