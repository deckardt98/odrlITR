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

# Kernel specifications are kept deliberately lightweight. Character kernels
# cover the common cases, while a function (or a list containing `fun`) lets a
# caller provide K(x, y) without adding a mandatory SVM dependency.
.odrl_resolve_kernel <- function(kernel, control = NULL, candidate = list()) {
  if (.odrl_is_series_kernel(kernel)) {
    return(.odrl_resolve_series_kernel(kernel))
  }
  if (is.function(kernel)) {
    return(list(
      name = "custom", fun = kernel, args = list(), builtin = FALSE
    ))
  }
  if (is.list(kernel)) {
    listed_name <- kernel$name %||% NULL
    if (!is.null(listed_name) && length(listed_name) == 1L &&
        !is.na(listed_name)) {
      listed_name <- tolower(as.character(listed_name))
      listed_name <- switch(listed_name,
        gaussian = "rbf", poly = "polynomial", listed_name
      )
    }
    if (!is.null(listed_name) &&
        listed_name %in% c("linear", "rbf", "polynomial") &&
        is.null(kernel$fun) && is.null(kernel$kernel)) {
      if (listed_name != "polynomial") {
        return(list(name = listed_name, builtin = TRUE))
      }
      degree_grid <- kernel$degree %||%
        control$svm_polynomial_degree %||% 2L
      degree <- candidate$degree %||% degree_grid[[1L]]
      scale <- kernel$scale %||% control$svm_polynomial_scale %||% 1
      offset <- kernel$offset %||% control$svm_polynomial_offset %||% 1
      if (length(degree) != 1L || !is.finite(degree) || degree < 1 ||
          degree != as.integer(degree)) {
        .odrl_abort("Polynomial kernel degree must be one positive integer.")
      }
      .odrl_check_scalar(scale, "svm_polynomial_scale", 0, Inf,
                         open_lower = TRUE)
      .odrl_check_scalar(offset, "svm_polynomial_offset", 0, Inf)
      return(list(
        name = "polynomial", degree = as.integer(degree), scale = scale,
        offset = offset, builtin = TRUE
      ))
    }
    fun <- kernel$fun %||% kernel$kernel
    if (!is.function(fun)) {
      .odrl_abort(
        "A custom `svm_kernel` list must contain a function named `fun` ",
        "(or `kernel`) that returns K(x, y)."
      )
    }
    args <- kernel$args %||% list()
    if (!is.list(args)) .odrl_abort("Custom kernel `args` must be a list.")
    name <- as.character(kernel$name %||% "custom")
    if (length(name) != 1L || is.na(name) || !nzchar(name)) {
      .odrl_abort("A custom kernel `name` must be one nonempty string.")
    }
    return(list(
      name = name,
      fun = fun,
      args = args,
      builtin = FALSE
    ))
  }
  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel)) {
    .odrl_abort(
      "`svm_kernel` must be a supported name, a kernel function, or a list ",
      "containing a kernel function."
    )
  }
  name <- tolower(kernel)
  name <- switch(name, gaussian = "rbf", poly = "polynomial", name)
  if (!name %in% c("linear", "rbf", "polynomial", "custom")) {
    .odrl_abort(
      "Unknown `svm_kernel`: ", kernel, ". Use `\"linear\"`, `\"rbf\"` ",
      "(`\"gaussian\"`), `\"polynomial\"`, or a custom kernel function."
    )
  }
  if (name == "custom") {
    fun <- control$svm_kernel_function %||% NULL
    if (!is.function(fun)) {
      .odrl_abort(
        "`svm_kernel = \"custom\"` requires `svm_kernel_function`."
      )
    }
    args <- control$svm_kernel_args %||% list()
    if (!is.list(args)) .odrl_abort("`svm_kernel_args` must be a list.")
    return(list(name = "custom", fun = fun, args = args, builtin = FALSE))
  }
  if (name == "polynomial") {
    degree_grid <- control$svm_polynomial_degree %||% 2L
    degree <- candidate$degree %||% degree_grid[[1L]]
    scale <- control$svm_polynomial_scale %||% 1
    offset <- control$svm_polynomial_offset %||% 1
    if (length(degree) != 1L || !is.finite(degree) || degree < 1 ||
        degree != as.integer(degree)) {
      .odrl_abort("Polynomial kernel degree must be one positive integer.")
    }
    .odrl_check_scalar(scale, "svm_polynomial_scale", 0, Inf,
                       open_lower = TRUE)
    .odrl_check_scalar(offset, "svm_polynomial_offset", 0, Inf)
    return(list(
      name = name, degree = as.integer(degree), scale = scale,
      offset = offset, builtin = TRUE
    ))
  }
  list(name = name, builtin = TRUE)
}

.odrl_builtin_kernel <- function(spec, name = NULL) {
  builtin <- isTRUE(spec$builtin) && is.null(spec$fun) &&
    !isTRUE(spec$finite_features)
  if (is.null(name)) builtin else builtin && identical(spec$name, name)
}

.odrl_effective_hinge_mode <- function(spec, mode = "auto") {
  if (identical(mode, "auto")) {
    if (.odrl_builtin_kernel(spec, "rbf")) "bounded" else "regularized"
  } else {
    mode
  }
}

.odrl_call_custom_kernel <- function(spec, x, y) {
  value <- tryCatch(
    do.call(spec$fun, c(list(x = x, y = y), spec$args)),
    error = function(error) {
      .odrl_abort("The custom kernel failed: ", conditionMessage(error))
    }
  )
  value <- as.matrix(value)
  if (!identical(dim(value), c(nrow(x), nrow(y))) ||
      !is.numeric(value) || any(!is.finite(value))) {
    .odrl_abort(
      "A custom kernel must return a finite numeric matrix with ",
      "`nrow(x)` rows and `nrow(y)` columns."
    )
  }
  value
}

.odrl_kernel_matrix <- function(x, y = x, kernel, bandwidth2 = NULL,
                                control = NULL, candidate = list()) {
  spec <- .odrl_resolve_kernel(
    kernel, control = control, candidate = candidate
  )
  if (inherits(spec, "odrl_series_kernel")) {
    .odrl_abort(
      "Finite-series specifications use a primal feature-map fit and do not ",
      "form a dense kernel matrix."
    )
  }
  if (!.odrl_builtin_kernel(spec)) {
    return(.odrl_call_custom_kernel(spec, x, y))
  }
  if (.odrl_builtin_kernel(spec, "linear")) {
    return(tcrossprod(x, y) / max(1, ncol(x)))
  }
  if (.odrl_builtin_kernel(spec, "rbf")) {
    if (length(bandwidth2) != 1L || !is.finite(bandwidth2) ||
        bandwidth2 <= 0) {
      .odrl_abort("A positive finite RBF squared bandwidth is required.")
    }
    return(exp(-.odrl_squared_distance(x, y) / (2 * bandwidth2)))
  }
  if (.odrl_builtin_kernel(spec, "polynomial")) {
    inner <- tcrossprod(x, y) / max(1, ncol(x))
    return((spec$scale * inner + spec$offset)^spec$degree)
  }
  .odrl_abort("Unsupported built-in kernel specification.")
}

.odrl_resolve_surrogate_loss <- function(loss) {
  if (inherits(loss, "odrl_surrogate_loss")) {
    return(loss)
  }
  if (is.character(loss) && length(loss) == 1L && !is.na(loss)) {
    name <- tolower(gsub("-", "_", loss, fixed = TRUE))
    name <- switch(name,
      logit = "logistic",
      squaredhinge = "squared_hinge",
      sq_hinge = "squared_hinge",
      name
    )
    if (name == "logistic") {
      return(structure(list(
        name = name,
        value = function(margin) .odrl_log1pexp(-margin),
        gradient = function(margin) -stats::plogis(-margin),
        builtin = TRUE
      ), class = c("odrl_surrogate_loss", "list")))
    }
    if (name == "exponential") {
      return(structure(list(
        name = name,
        value = function(margin) .odrl_exponential_margin(margin)$loss,
        gradient = function(margin) {
          .odrl_exponential_margin(margin)$gradient
        },
        builtin = TRUE
      ), class = c("odrl_surrogate_loss", "list")))
    }
    if (name == "hinge") {
      return(structure(list(
        name = name,
        value = function(margin) pmax(1 - margin, 0),
        gradient = function(margin) ifelse(margin < 1, -1, 0),
        builtin = TRUE
      ), class = c("odrl_surrogate_loss", "list")))
    }
    if (name == "squared_hinge") {
      return(structure(list(
        name = name,
        value = function(margin) pmax(1 - margin, 0)^2,
        gradient = function(margin) -2 * pmax(1 - margin, 0),
        builtin = TRUE
      ), class = c("odrl_surrogate_loss", "list")))
    }
    .odrl_abort(
      "Unknown SVM surrogate loss: ", loss, ". Use `\"hinge\"`, ",
      "`\"exponential\"`, `\"logistic\"`, `\"squared_hinge\"`, or a ",
      "custom loss specification."
    )
  }
  if (is.function(loss)) {
    gradient <- attr(loss, "gradient", exact = TRUE)
    if (!is.function(gradient)) {
      .odrl_abort(
        "A custom loss function must have a function-valued `gradient` ",
        "attribute giving the derivative (or a subgradient) with respect ",
        "to the signed margin."
      )
    }
    loss <- list(
      name = attr(loss, "name", exact = TRUE) %||% "custom",
      value = loss,
      gradient = gradient
    )
  }
  if (!is.list(loss)) {
    .odrl_abort(
      "A custom SVM loss must be a list with `value` and `gradient` functions."
    )
  }
  value <- loss$value %||% loss$loss
  gradient <- loss$gradient %||% loss$subgradient
  if (!is.function(value) || !is.function(gradient)) {
    .odrl_abort(
      "A custom SVM loss must contain function-valued `value` and ",
      "`gradient` (or `subgradient`) entries. Both operate on signed margins."
    )
  }
  probe <- c(-1, 0, 1)
  probe_value <- value(probe)
  probe_gradient <- gradient(probe)
  if (!is.numeric(probe_value) || length(probe_value) != length(probe) ||
      any(!is.finite(probe_value)) || any(probe_value < 0)) {
    .odrl_abort(
      "The custom loss `value` function must return one finite nonnegative ",
      "number per signed margin."
    )
  }
  if (!is.numeric(probe_gradient) ||
      length(probe_gradient) != length(probe) ||
      any(!is.finite(probe_gradient))) {
    .odrl_abort(
      "The custom loss `gradient` function must return one finite number ",
      "per signed margin."
    )
  }
  name <- as.character(loss$name %||% "custom")
  if (length(name) != 1L || is.na(name) || !nzchar(name)) {
    .odrl_abort("A custom loss `name` must be one nonempty string.")
  }
  structure(list(
    name = name,
    value = value,
    gradient = gradient,
    builtin = FALSE
  ), class = c("odrl_surrogate_loss", "list"))
}

.odrl_finalize_bounded_hinge <- function(fit) {
  raw <- fit$fitted
  fit$fitted_unclipped <- raw
  fit$fitted <- .odrl_hardtanh(raw)
  fit$globally_bounded <- TRUE
  fit$bounded_output <- TRUE
  fit$hinge_mode <- "bounded"
  fit$clipping <- "hard_tanh"
  fit
}

.odrl_fit_kernel_hinge <- function(x, score, kernel, multiplier, lambda,
                                   maxit, seed, control = NULL,
                                   candidate = list()) {
  spec <- .odrl_resolve_kernel(
    kernel, control = control, candidate = candidate
  )
  bandwidth2 <- if (.odrl_builtin_kernel(spec, "rbf")) {
    .odrl_median_squared_distance(x, seed) * multiplier
  } else {
    NA_real_
  }
  k <- .odrl_kernel_matrix(x, kernel = spec, bandwidth2 = bandwidth2)
  if (nrow(k) == ncol(k) &&
      max(abs(k - t(k))) > 1e-7 * max(1, max(abs(k)))) {
    .odrl_abort("The training kernel matrix must be symmetric.")
  }
  label <- ifelse(score >= 0, 1, -1)
  scale <- mean(abs(score))
  weight <- if (is.finite(scale) && scale > 0) {
    abs(score) / scale
  } else {
    rep(1, length(score))
  }
  upper <- weight / length(score)
  active <- which(upper > 0)
  # Penalizing the constant offset is equivalent to using the kernel K + 1.
  augmented <- k + 1
  dual_fit <- NULL
  u <- numeric(length(score))
  if (length(active)) {
    active_kernel <- augmented[active, active, drop = FALSE]
    active_label <- label[active]
    objective <- function(active_u) {
      signed <- active_u * active_label
      0.5 * drop(crossprod(signed, active_kernel %*% signed)) / lambda -
        sum(active_u)
    }
    gradient <- function(active_u) {
      signed <- active_u * active_label
      active_label * drop(active_kernel %*% signed) / lambda - 1
    }
    set.seed(seed)
    dual_fit <- stats::optim(
      par = upper[active] / 2, fn = objective, gr = gradient,
      method = "L-BFGS-B", lower = 0, upper = upper[active],
      control = list(maxit = maxit, factr = 1e8)
    )
    u[active] <- pmin(pmax(dual_fit$par, 0), upper[active])
  }
  signed <- u * label
  intercept <- sum(signed) / lambda
  alpha <- signed / lambda
  fitted <- intercept + drop(k %*% alpha)
  norm_squared <- max(
    intercept^2 + drop(crossprod(alpha, k %*% alpha)), 0
  )
  primal <- mean(weight * pmax(1 - label * fitted, 0)) +
    0.5 * lambda * norm_squared
  dual <- sum(u) -
    0.5 * drop(crossprod(signed, augmented %*% signed)) / lambda
  fit <- list(
    intercept = intercept,
    alpha = alpha,
    training_x = x,
    kernel_spec = spec,
    bandwidth2 = bandwidth2,
    lambda = lambda,
    loss = "hinge",
    loss_spec = .odrl_resolve_surrogate_loss("hinge"),
    fitted = fitted,
    convergence = dual_fit$convergence %||% 0L,
    message = dual_fit$message %||% NULL,
    attempts = 1L,
    objective = primal,
    dual_objective = dual,
    duality_gap = max(primal - dual, 0),
    rkhs_norm = sqrt(norm_squared),
    globally_bounded = FALSE,
    bounded_output = FALSE,
    hinge_mode = "regularized"
  )
  fit
}

.odrl_fit_kernel_surrogate <- function(x, score, kernel, multiplier, lambda,
                                       maxit, seed, loss, control = NULL,
                                       candidate = list()) {
  spec <- .odrl_resolve_kernel(kernel, control = control,
                               candidate = candidate)
  bandwidth2 <- if (.odrl_builtin_kernel(spec, "rbf")) {
    .odrl_median_squared_distance(x, seed) * multiplier
  } else {
    NA_real_
  }
  k <- .odrl_kernel_matrix(x, kernel = spec, bandwidth2 = bandwidth2)
  if (nrow(k) == ncol(k) &&
      max(abs(k - t(k))) > 1e-7 * max(1, max(abs(k)))) {
    .odrl_abort("The training kernel matrix must be symmetric.")
  }
  loss_spec <- .odrl_resolve_surrogate_loss(loss)
  label <- ifelse(score >= 0, 1, -1)
  score_scale <- mean(abs(score))
  weight <- if (is.finite(score_scale) && score_scale > 0) {
    abs(score) / score_scale
  } else {
    rep(1, length(score))
  }
  n <- length(score)
  objective <- function(theta) {
    intercept <- theta[[1L]]
    alpha <- theta[-1L]
    raw <- intercept + drop(k %*% alpha)
    margin <- label * raw
    mean(weight * loss_spec$value(margin)) +
      0.5 * lambda * (intercept^2 +
        drop(crossprod(alpha, k %*% alpha)))
  }
  gradient <- function(theta) {
    intercept <- theta[[1L]]
    alpha <- theta[-1L]
    raw <- intercept + drop(k %*% alpha)
    margin <- label * raw
    derivative <- weight * loss_spec$gradient(margin) * label / n
    c(sum(derivative) + lambda * intercept,
      drop(k %*% derivative) + lambda * drop(k %*% alpha))
  }
  set.seed(seed)
  initial_intercept <- if (isTRUE(loss_spec$builtin) &&
      identical(loss_spec$name, "exponential")) {
    positive <- max(sum(weight[label == 1]), .Machine$double.eps)
    negative <- max(sum(weight[label == -1]), .Machine$double.eps)
    0.5 * log(positive / negative)
  } else {
    stats::qlogis(pmin(pmax(mean(label == 1), 0.01), 0.99))
  }
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
    kernel_spec = spec,
    bandwidth2 = bandwidth2,
    lambda = lambda,
    loss = loss_spec$name,
    loss_spec = loss_spec,
    fitted = fit$par[[1L]] + drop(k %*% fit$par[-1L]),
    convergence = fit$convergence,
    message = fit$message,
    attempts = attempts,
    objective = fit$value,
    rkhs_norm = sqrt(max(
      fit$par[[1L]]^2 +
        drop(crossprod(fit$par[-1L], k %*% fit$par[-1L])),
      0
    )),
    globally_bounded = FALSE,
    bounded_output = FALSE,
    hinge_mode = if (isTRUE(loss_spec$builtin) &&
      identical(loss_spec$name, "hinge")) "regularized" else NULL
  )
}

.odrl_predict_kernel_unclipped <- function(fit, newx, kernel = NULL) {
  if (!is.null(fit$series_map)) {
    return(.odrl_predict_series_unclipped(fit, newx))
  }
  spec <- fit$kernel_spec %||% .odrl_resolve_kernel(kernel)
  k <- .odrl_kernel_matrix(
    newx, fit$training_x, kernel = spec, bandwidth2 = fit$bandwidth2
  )
  (fit$intercept %||% 0) + drop(k %*% fit$alpha)
}

.odrl_predict_kernel_raw <- function(fit, newx, kernel = NULL) {
  raw <- .odrl_predict_kernel_unclipped(fit, newx, kernel)
  if (isTRUE(fit$bounded_output)) {
    .odrl_hardtanh(raw)
  } else {
    raw
  }
}

.odrl_kernel_grid <- function(control) {
  kernel_spec <- .odrl_resolve_kernel(control$svm_kernel, control = control)
  multiplier <- if (.odrl_builtin_kernel(kernel_spec, "rbf")) {
    control$svm_rbf_multiplier
  } else {
    1
  }
  degree <- if (.odrl_builtin_kernel(kernel_spec, "polynomial")) {
    control$svm_polynomial_degree %||% 2L
  } else {
    NA_integer_
  }
  base <- if (inherits(kernel_spec, "odrl_series_kernel")) {
    .odrl_series_grid(kernel_spec)
  } else if (.odrl_builtin_kernel(kernel_spec, "polynomial")) {
    expand.grid(
      multiplier = multiplier,
      degree = degree,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(multiplier = multiplier)
  }
  merge(
    base, data.frame(lambda = rev(control$svm_penalty)), by = NULL,
    sort = FALSE
  )
}

.odrl_fit_svm_candidate <- function(x, score, control, loss, candidate, seed) {
  loss_spec <- .odrl_resolve_surrogate_loss(loss)
  kernel_spec <- .odrl_resolve_kernel(
    control$svm_kernel, control = control, candidate = candidate
  )
  hinge_mode <- .odrl_effective_hinge_mode(
    kernel_spec, control$svm_hinge_mode %||% "auto"
  )
  bounded_hinge <- isTRUE(loss_spec$builtin) &&
    identical(loss_spec$name, "hinge") &&
    identical(hinge_mode, "bounded")
  if (inherits(kernel_spec, "odrl_series_kernel")) {
    fit <- .odrl_fit_series_surrogate(
      x = x, score = score, spec = kernel_spec, candidate = candidate,
      lambda = candidate$lambda, maxit = control$svm_maxit, seed = seed,
      loss = loss_spec
    )
    if (bounded_hinge) .odrl_finalize_bounded_hinge(fit) else fit
  } else if (isTRUE(loss_spec$builtin) &&
      identical(loss_spec$name, "hinge")) {
    fit <- .odrl_fit_kernel_hinge(
      x, score, kernel_spec, candidate$multiplier,
      candidate$lambda, control$svm_maxit, seed,
      control = control, candidate = candidate
    )
    if (bounded_hinge) .odrl_finalize_bounded_hinge(fit) else fit
  } else {
    .odrl_fit_kernel_surrogate(
      x, score, kernel_spec, candidate$multiplier,
      candidate$lambda, control$svm_maxit, seed, loss = loss_spec,
      control = control, candidate = candidate
    )
  }
}

#' Fit a tuned kernel ODRL rule
#' @noRd
.odrl_fit_svm <- function(x, score, control, loss) {
  loss_spec <- .odrl_resolve_surrogate_loss(loss)
  kernel_spec <- .odrl_resolve_kernel(control$svm_kernel, control = control)
  hinge_mode <- .odrl_effective_hinge_mode(
    kernel_spec, control$svm_hinge_mode %||% "auto"
  )
  folds <- min(control$svm_folds, nrow(x))
  fold_id <- .odrl_score_folds(score, folds, control$seed + 300L)
  grid <- .odrl_kernel_grid(control)
  criterion <- matrix(NA_real_, nrow(grid), folds)
  convergence <- matrix(NA_integer_, nrow(grid), folds)
  started <- proc.time()[["elapsed"]]
  for (j in seq_len(nrow(grid))) {
    candidate <- as.list(grid[j, , drop = FALSE])
    for (fold in seq_len(folds)) {
      train <- fold_id != fold
      holdout <- !train
      if (inherits(kernel_spec, "odrl_series_kernel")) {
        train_x <- x[train, , drop = FALSE]
        holdout_x <- x[holdout, , drop = FALSE]
      } else {
        transform <- .odrl_standardize_fit(x[train, , drop = FALSE])
        train_x <- .odrl_standardize_apply(
          x[train, , drop = FALSE], transform
        )
        holdout_x <- .odrl_standardize_apply(
          x[holdout, , drop = FALSE], transform
        )
      }
      fit <- .odrl_fit_svm_candidate(
        train_x, score[train], control, loss_spec, candidate,
        control$seed + 10000L * j + fold
      )
      raw <- .odrl_predict_kernel_raw(fit, holdout_x)
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
  if (inherits(kernel_spec, "odrl_series_kernel")) {
    transform <- NULL
    standardized_x <- x
  } else {
    transform <- .odrl_standardize_fit(x)
    standardized_x <- .odrl_standardize_apply(x, transform)
  }
  final <- .odrl_fit_svm_candidate(
    standardized_x, score, control, loss_spec,
    as.list(grid[best, , drop = FALSE]), control$seed + 900000L
  )
  if (final$convergence != 0L) {
    .odrl_abort(
      "The selected kernel candidate did not converge on the full sample. ",
      "Increase `svm_maxit`."
    )
  }
  raw <- .odrl_predict_kernel_raw(final, standardized_x)
  action <- ifelse(raw >= 0, 1, -1)
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    engine = if (inherits(kernel_spec, "odrl_series_kernel")) {
      "series"
    } else {
      "kernel"
    },
    kernel = kernel_spec$name,
    kernel_spec = final$kernel_spec,
    loss = loss_spec$name,
    loss_spec = loss_spec,
    fit = final,
    transform = transform,
    tuning = grid,
    selected = grid[best, , drop = FALSE],
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      folds = folds,
      elapsed = elapsed,
      kernel = kernel_spec$name,
      loss = loss_spec$name,
      custom_kernel = !.odrl_builtin_kernel(kernel_spec) &&
        !inherits(kernel_spec, "odrl_series_kernel"),
      finite_series = inherits(kernel_spec, "odrl_series_kernel"),
      series_combine = if (inherits(kernel_spec, "odrl_series_kernel")) {
        kernel_spec$combine
      } else {
        NA_character_
      },
      series_features = final$series_map$feature_count %||% NA_integer_,
      custom_loss = !isTRUE(loss_spec$builtin),
      hinge_mode = final$hinge_mode %||% NA_character_,
      clipping = final$clipping %||% NA_character_,
      globally_bounded = isTRUE(final$globally_bounded),
      rkhs_norm = final$rkhs_norm %||% NA_real_,
      duality_gap = final$duality_gap %||% NA_real_,
      max_training_abs_unclipped = if (isTRUE(final$bounded_output)) {
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
  z <- if (is.null(object$transform)) {
    newx
  } else {
    .odrl_standardize_apply(newx, object$transform)
  }
  raw <- .odrl_predict_kernel_raw(object$fit, z)
  if (type == "score") raw else ifelse(raw >= 0, 1, -1)
}
