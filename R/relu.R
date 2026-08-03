.odrl_relu_unpack <- function(theta, p, hidden) {
  if (hidden == 0L) {
    return(list(intercept = theta[[1L]], beta = theta[-1L]))
  }
  end_w1 <- p * hidden
  end_b1 <- end_w1 + hidden
  end_w2 <- end_b1 + hidden
  list(
    w1 = matrix(theta[seq_len(end_w1)], nrow = p, ncol = hidden),
    b1 = theta[(end_w1 + 1L):end_b1],
    w2 = theta[(end_b1 + 1L):end_w2],
    b2 = theta[[end_w2 + 1L]]
  )
}

.odrl_relu_forward <- function(theta, x, hidden) {
  p <- ncol(x)
  parameter <- .odrl_relu_unpack(theta, p, hidden)
  if (hidden == 0L) {
    raw <- parameter$intercept + drop(x %*% parameter$beta)
    return(list(raw = raw, parameter = parameter))
  }
  pre <- sweep(x %*% parameter$w1, 2L, parameter$b1, "+")
  activation <- pmax(pre, 0)
  raw <- parameter$b2 + drop(activation %*% parameter$w2)
  list(
    raw = raw, pre = pre, activation = activation, parameter = parameter
  )
}

.odrl_relu_objective <- function(theta, x, label, weight, hidden, decay,
                                 loss, gradient = FALSE) {
  n <- nrow(x)
  forward <- .odrl_relu_forward(theta, x, hidden)
  raw <- forward$raw
  if (loss == "hinge") {
    bounded <- .odrl_hardtanh(raw)
    value <- mean(weight * (1 - label * bounded))
    derivative <- -weight * label * as.numeric(abs(raw) < 1) / n
  } else {
    value <- mean(weight * .odrl_log1pexp(-label * raw))
    derivative <- -weight * label * stats::plogis(-label * raw) / n
  }
  parameter <- forward$parameter
  if (hidden == 0L) {
    penalty <- 0.5 * decay * sum(parameter$beta^2)
    if (!gradient) return(value + penalty)
    return(c(sum(derivative), drop(crossprod(x, derivative)) +
               decay * parameter$beta))
  }
  penalty <- 0.5 * decay *
    (sum(parameter$w1^2) + sum(parameter$w2^2))
  if (!gradient) return(value + penalty)
  hidden_derivative <- tcrossprod(derivative, parameter$w2) *
    (forward$pre > 0)
  grad_w1 <- crossprod(x, hidden_derivative) + decay * parameter$w1
  grad_b1 <- colSums(hidden_derivative)
  grad_w2 <- drop(crossprod(forward$activation, derivative)) +
    decay * parameter$w2
  grad_b2 <- sum(derivative)
  c(as.vector(grad_w1), grad_b1, grad_w2, grad_b2)
}

.odrl_relu_initial <- function(x, label, hidden, seed, loss) {
  p <- ncol(x)
  intercept <- if (loss == "hinge") {
    0
  } else {
    stats::qlogis(pmin(pmax(mean(label == 1), 0.01), 0.99))
  }
  if (hidden == 0L) return(c(intercept, rep(0, p)))
  set.seed(seed)
  scale <- sqrt(2 / max(1, p))
  c(
    stats::rnorm(p * hidden, sd = scale),
    rep(0, hidden),
    stats::rnorm(hidden, sd = 0.05),
    intercept
  )
}

.odrl_fit_relu_candidate <- function(x, score, hidden, decay, loss,
                                     restarts, maxit, seed) {
  label <- ifelse(score >= 0, 1, -1)
  weight <- abs(score) / mean(abs(score))
  fits <- vector("list", restarts)
  for (restart in seq_len(restarts)) {
    initial <- .odrl_relu_initial(
      x, label, hidden, seed + restart * 1009L, loss
    )
    fits[[restart]] <- stats::optim(
      par = initial,
      fn = function(theta) .odrl_relu_objective(
        theta, x, label, weight, hidden, decay, loss, gradient = FALSE
      ),
      gr = function(theta) .odrl_relu_objective(
        theta, x, label, weight, hidden, decay, loss, gradient = TRUE
      ),
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e8)
    )
  }
  values <- vapply(fits, `[[`, numeric(1), "value")
  converged <- which(vapply(fits, `[[`, integer(1), "convergence") == 0L)
  if (!length(converged)) {
    provisional <- which.min(values)
    retry <- stats::optim(
      par = fits[[provisional]]$par,
      fn = function(theta) .odrl_relu_objective(
        theta, x, label, weight, hidden, decay, loss, gradient = FALSE
      ),
      gr = function(theta) .odrl_relu_objective(
        theta, x, label, weight, hidden, decay, loss, gradient = TRUE
      ),
      method = "L-BFGS-B",
      control = list(maxit = max(3L * maxit, maxit + 100L), factr = 1e8)
    )
    fits[[length(fits) + 1L]] <- retry
    values <- c(values, retry$value)
    if (retry$convergence == 0L) converged <- length(fits)
  }
  best <- if (length(converged)) {
    converged[[which.min(values[converged])]]
  } else {
    which.min(values)
  }
  fit <- fits[[best]]
  list(
    theta = fit$par,
    hidden = hidden,
    decay = decay,
    loss = loss,
    objective = fit$value,
    convergence = fit$convergence,
    message = fit$message,
    restart = best,
    attempts = length(fits),
    all_objectives = values
  )
}

.odrl_predict_relu_raw <- function(fit, newx) {
  raw <- .odrl_relu_forward(fit$theta, newx, fit$hidden)$raw
  if (fit$loss == "hinge") .odrl_hardtanh(raw) else raw
}

.odrl_select_relu <- function(grid, selection) {
  converged <- which(
    grid$convergence_failures == 0L & is.finite(grid$mean_criterion)
  )
  if (!length(converged)) {
    .odrl_abort(
      "Every ReLU candidate had at least one nonconverged tuning fit. ",
      "Increase `relu_maxit`, add restarts, or simplify the tuning grid."
    )
  }
  best <- converged[[which.max(grid$mean_criterion[converged])]]
  if (selection == "best") return(best)
  cutoff <- grid$mean_criterion[[best]] - grid$se_criterion[[best]]
  eligible <- converged[grid$mean_criterion[converged] >= cutoff]
  ordered <- eligible[order(
    grid$hidden[eligible], -grid$decay[eligible],
    -grid$mean_criterion[eligible]
  )]
  ordered[[1L]]
}

#' Fit a tuned ReLU ODRL rule
#' @noRd
.odrl_fit_relu <- function(x, score, control, loss) {
  folds <- min(control$relu_folds, nrow(x))
  fold_id <- .odrl_score_folds(score, folds, control$seed + 400L)
  grid <- expand.grid(
    hidden = control$relu_hidden_units,
    decay = rev(control$relu_decay),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  criterion <- matrix(NA_real_, nrow(grid), folds)
  convergence <- matrix(NA_integer_, nrow(grid), folds)
  started <- proc.time()[["elapsed"]]
  for (j in seq_len(nrow(grid))) {
    for (fold in seq_len(folds)) {
      train <- fold_id != fold
      holdout <- !train
      transform <- .odrl_minmax_fit(x[train, , drop = FALSE])
      train_x <- .odrl_minmax_apply(x[train, , drop = FALSE], transform)
      holdout_x <- .odrl_minmax_apply(x[holdout, , drop = FALSE], transform)
      fit <- .odrl_fit_relu_candidate(
        train_x, score[train], grid$hidden[[j]], grid$decay[[j]], loss,
        control$relu_restarts, control$relu_maxit,
        control$seed + 20000L * j + fold
      )
      raw <- .odrl_predict_relu_raw(fit, holdout_x)
      action <- ifelse(raw >= 0, 1, -1)
      criterion[j, fold] <- .odrl_empirical_criterion(
        score[holdout], action
      )
      convergence[j, fold] <- fit$convergence
    }
  }
  grid$mean_criterion <- rowMeans(criterion)
  grid$se_criterion <- apply(criterion, 1L, stats::sd) / sqrt(folds)
  grid$convergence_failures <- rowSums(convergence != 0L, na.rm = TRUE)
  selected <- .odrl_select_relu(grid, control$relu_selection)
  transform <- .odrl_minmax_fit(x)
  transformed_x <- .odrl_minmax_apply(x, transform)
  final <- .odrl_fit_relu_candidate(
    transformed_x, score, grid$hidden[[selected]], grid$decay[[selected]],
    loss, control$relu_refit_restarts, control$relu_maxit,
    control$seed + 950000L
  )
  if (final$convergence != 0L) {
    .odrl_abort(
      "The selected ReLU candidate did not converge on the full sample. ",
      "Increase `relu_maxit` or `relu_refit_restarts`."
    )
  }
  raw <- .odrl_predict_relu_raw(final, transformed_x)
  action <- ifelse(raw >= 0, 1, -1)
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    engine = "relu",
    loss = loss,
    fit = final,
    transform = transform,
    tuning = grid,
    selected = grid[selected, , drop = FALSE],
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      folds = folds,
      selection_rule = control$relu_selection,
      elapsed = elapsed,
      globally_bounded = identical(loss, "hinge"),
      optimization_score_scale = "mean-absolute-score normalized",
      convergence = final$convergence,
      selected_restart = final$restart
    )
  ), class = "odrl_policy_relu")
}

.odrl_predict_relu <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
  z <- .odrl_minmax_apply(newx, object$transform)
  raw <- .odrl_predict_relu_raw(object$fit, z)
  if (type == "score") raw else ifelse(raw >= 0, 1, -1)
}
