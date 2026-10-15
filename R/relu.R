.odrl_relu_architecture <- function(hidden) {
  hidden <- as.integer(hidden)
  if (!length(hidden) || (length(hidden) == 1L && hidden == 0L)) {
    return(integer())
  }
  if (anyNA(hidden) || any(hidden <= 0L)) {
    .odrl_abort(
      "A neural network architecture must be a vector of positive layer ",
      "widths; use `integer()` or `0L` for an affine score."
    )
  }
  hidden
}

.odrl_relu_architecture_label <- function(hidden) {
  hidden <- .odrl_relu_architecture(hidden)
  if (!length(hidden)) "affine" else paste(hidden, collapse = "x")
}

.odrl_relu_activation <- function(x, activation, leaky_slope) {
  switch(
    activation,
    relu = pmax(x, 0),
    leaky_relu = ifelse(x >= 0, x, leaky_slope * x),
    tanh = tanh(x),
    sigmoid = stats::plogis(x),
    linear = x,
    .odrl_abort("Unsupported neural network activation `", activation, "`.")
  )
}

.odrl_relu_activation_derivative <- function(x, activation, leaky_slope) {
  switch(
    activation,
    relu = as.numeric(x > 0),
    leaky_relu = ifelse(x >= 0, 1, leaky_slope),
    tanh = 1 - tanh(x)^2,
    sigmoid = {
      probability <- stats::plogis(x)
      probability * (1 - probability)
    },
    linear = rep(1, length(x)),
    .odrl_abort("Unsupported neural network activation `", activation, "`.")
  )
}

.odrl_relu_unpack <- function(theta, p, hidden) {
  architecture <- .odrl_relu_architecture(hidden)
  if (!length(architecture)) {
    return(list(intercept = theta[[1L]], beta = theta[-1L]))
  }
  dimensions <- c(p, architecture, 1L)
  weights <- biases <- vector("list", length(dimensions) - 1L)
  cursor <- 0L
  # Each layer contributes its weight matrix, stored column by column,
  # followed by its bias vector. Initialization and gradients use this order.
  for (layer in seq_along(weights)) {
    weight_count <- dimensions[[layer]] * dimensions[[layer + 1L]]
    weight_index <- cursor + seq_len(weight_count)
    weights[[layer]] <- matrix(
      theta[weight_index], nrow = dimensions[[layer]],
      ncol = dimensions[[layer + 1L]]
    )
    cursor <- cursor + weight_count
    bias_index <- cursor + seq_len(dimensions[[layer + 1L]])
    biases[[layer]] <- theta[bias_index]
    cursor <- cursor + dimensions[[layer + 1L]]
  }
  list(weights = weights, biases = biases, architecture = architecture)
}

.odrl_relu_forward <- function(theta, x, hidden, activation = "relu",
                               leaky_slope = 0.01) {
  p <- ncol(x)
  parameter <- .odrl_relu_unpack(theta, p, hidden)
  architecture <- .odrl_relu_architecture(hidden)
  if (!length(architecture)) {
    raw <- parameter$intercept + drop(x %*% parameter$beta)
    return(list(raw = raw, parameter = parameter, inputs = list(x)))
  }
  # Save each layer's input and its value before activation for backpropagation.
  inputs <- list(x)
  pre_activation <- activations <- vector("list", length(architecture))
  layer_input <- x
  for (layer in seq_along(architecture)) {
    pre_activation[[layer]] <- sweep(
      layer_input %*% parameter$weights[[layer]], 2L,
      parameter$biases[[layer]], "+"
    )
    layer_output <- matrix(
      .odrl_relu_activation(pre_activation[[layer]], activation, leaky_slope),
      nrow = nrow(x), ncol = architecture[[layer]]
    )
    activations[[layer]] <- layer_output
    inputs[[layer + 1L]] <- layer_output
    layer_input <- layer_output
  }
  # The output layer is affine. Hinge clipping is applied by the loss and
  # prediction functions so that the unclipped score remains available here.
  output_layer <- length(parameter$weights)
  raw <- parameter$biases[[output_layer]] +
    drop(layer_input %*% parameter$weights[[output_layer]])
  list(
    raw = raw, pre = pre_activation, activations = activations, inputs = inputs,
    parameter = parameter
  )
}

.odrl_relu_custom_loss <- function(loss, margin) {
  result <- loss(margin)
  if (!is.list(result) ||
      !all(c("loss", "gradient") %in% names(result))) {
    .odrl_abort(
      "A custom surrogate loss must return a list with numeric `loss` and ",
      "`gradient` vectors when called with the signed margin."
    )
  }
  value <- result$loss
  gradient <- result$gradient
  if (!is.numeric(value) || !is.numeric(gradient) ||
      length(value) != length(margin) ||
      length(gradient) != length(margin) ||
      any(!is.finite(value)) || any(value < 0) ||
      any(!is.finite(gradient))) {
    .odrl_abort(
      "A custom surrogate loss returned invalid values. `loss` and ",
      "`gradient` must be finite numeric vectors with one value per margin, ",
      "and `loss` must be nonnegative."
    )
  }
  list(loss = value, gradient = gradient)
}

.odrl_relu_loss_name <- function(loss) {
  if (is.function(loss)) {
    name <- attr(loss, "name", exact = TRUE)
    if (is.character(name) && length(name) == 1L && nzchar(name)) name else
      "custom"
  } else {
    as.character(loss)
  }
}

.odrl_relu_objective <- function(theta, x, label, weight, hidden, decay,
                                 loss, gradient = FALSE,
                                 activation = "relu",
                                 leaky_slope = 0.01) {
  n <- nrow(x)
  forward <- .odrl_relu_forward(
    theta, x, hidden, activation = activation, leaky_slope = leaky_slope
  )
  raw <- forward$raw
  # Minimize weighted mean loss plus a penalty on weights.
  # Biases are unpenalized.
  # score_derivative is the derivative of the mean loss with respect to each
  # unclipped network score, before the chain rule passes through the layers.
  if (identical(loss, "hinge")) {
    # Clipping to [-1, 1] makes hinge loss equal to 1 - label * bounded.
    # Outside the interval, clipping has zero derivative.
    bounded <- .odrl_hardtanh(raw)
    value <- mean(weight * (1 - label * bounded))
    score_derivative <- -weight * label * as.numeric(abs(raw) < 1) / n
  } else {
    # Positive margins mean agreement with the score's preferred treatment.
    margin <- label * raw
    evaluated <- if (identical(loss, "logistic")) {
      list(
        loss = .odrl_log1pexp(-margin),
        gradient = -stats::plogis(-margin)
      )
    } else if (identical(loss, "exponential")) {
      .odrl_exponential_margin(margin)
    } else if (identical(loss, "squared_hinge")) {
      violation <- pmax(0, 1 - margin)
      list(loss = violation^2, gradient = -2 * violation)
    } else if (is.function(loss)) {
      .odrl_relu_custom_loss(loss, margin)
    } else {
      .odrl_abort("Unsupported neural network surrogate loss.")
    }
    value <- mean(weight * evaluated$loss)
    score_derivative <- weight * evaluated$gradient * label / n
  }
  parameter <- forward$parameter
  architecture <- .odrl_relu_architecture(hidden)
  if (!length(architecture)) {
    penalty <- 0.5 * decay * sum(parameter$beta^2)
    if (!gradient) return(value + penalty)
    return(c(sum(score_derivative), drop(crossprod(x, score_derivative)) +
               decay * parameter$beta))
  }
  penalty <- 0.5 * decay * sum(vapply(
    parameter$weights, function(weight_matrix) sum(weight_matrix^2), numeric(1)
  ))
  if (!gradient) return(value + penalty)

  layer_count <- length(parameter$weights)
  weight_gradients <- bias_gradients <- vector("list", layer_count)
  layer_derivative <- matrix(score_derivative, ncol = 1L)
  # Walk backward from the output. At each layer, accumulate parameter
  # gradients, then propagate the derivative to the preceding hidden layer.
  for (layer in rev(seq_len(layer_count))) {
    weight_gradients[[layer]] <-
      crossprod(forward$inputs[[layer]], layer_derivative) +
        decay * parameter$weights[[layer]]
    bias_gradients[[layer]] <- colSums(layer_derivative)
    if (layer > 1L) {
      layer_derivative <- layer_derivative %*% t(parameter$weights[[layer]])
      activation_derivative <- matrix(
        .odrl_relu_activation_derivative(
          forward$pre[[layer - 1L]], activation, leaky_slope
        ),
        nrow = n, ncol = architecture[[layer - 1L]]
      )
      layer_derivative <- layer_derivative * activation_derivative
    }
  }
  unlist(Map(
    function(weight_gradient, bias_gradient) {
      c(as.vector(weight_gradient), bias_gradient)
    },
    weight_gradients, bias_gradients
  ), use.names = FALSE)
}

.odrl_relu_initial <- function(x, label, hidden, seed, loss,
                               activation = "relu", weight = NULL) {
  p <- ncol(x)
  architecture <- .odrl_relu_architecture(hidden)
  intercept <- if (identical(loss, "logistic")) {
    stats::qlogis(pmin(pmax(mean(label == 1), 0.01), 0.99))
  } else if (identical(loss, "exponential")) {
    if (is.null(weight)) weight <- rep(1, length(label))
    positive <- max(sum(weight[label == 1]), .Machine$double.eps)
    negative <- max(sum(weight[label == -1]), .Machine$double.eps)
    0.5 * log(positive / negative)
  } else {
    0
  }
  if (!length(architecture)) return(c(intercept, rep(0, p)))
  set.seed(seed)
  dimensions <- c(p, architecture, 1L)
  theta <- numeric()
  for (layer in seq_len(length(dimensions) - 1L)) {
    input_width <- dimensions[[layer]]
    output_width <- dimensions[[layer + 1L]]
    # Start the output near its intercept; scale hidden weights by layer width.
    weight_sd <- if (layer == length(dimensions) - 1L) {
      0.05
    } else if (activation %in% c("relu", "leaky_relu")) {
      sqrt(2 / max(1, input_width))
    } else {
      sqrt(1 / max(1, input_width))
    }
    theta <- c(
      theta,
      stats::rnorm(input_width * output_width, sd = weight_sd),
      if (layer == length(dimensions) - 1L) intercept else
        rep(0, output_width)
    )
  }
  theta
}

.odrl_relu_backend_spec <- function(backend, options = list()) {
  if (is.list(backend) && isTRUE(backend$native) &&
      identical(backend$name %||% "native", "native")) {
    if (length(options)) {
      .odrl_abort(
        "`relu_backend_options` are only used by a named external backend."
      )
    }
    return(backend)
  }
  if (is.null(backend) || identical(backend, "native")) {
    if (length(options)) {
      .odrl_abort(
        "`relu_backend_options` are only used by a named external backend."
      )
    }
    return(list(name = "native", native = TRUE))
  }
  if (identical(backend, "nnet")) {
    specification <- .odrl_nnet_backend(options)
    return(list(
      name = specification$name, native = FALSE,
      fit = specification$fit, predict = specification$predict
    ))
  }
  if (!is.list(backend) || !is.function(backend$fit) ||
      !is.function(backend$predict)) {
    .odrl_abort(
      "`relu_backend` must be `\"native\"`, `\"nnet\"`, or a custom ",
      "list containing `fit` and `predict` functions."
    )
  }
  if (length(options)) {
    .odrl_abort(
      "For a custom neural backend, close over options for the backend in ",
      "its `fit` and `predict` callbacks."
    )
  }
  name <- backend$name %||% "custom"
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    .odrl_abort("A custom neural backend's `name` must be one string.")
  }
  list(name = name, native = FALSE, fit = backend$fit,
       predict = backend$predict)
}

.odrl_fit_relu_candidate <- function(x, score, hidden, decay, loss,
                                     restarts, maxit, seed,
                                     activation = "relu",
                                     leaky_slope = 0.01,
                                     backend = "native") {
  architecture <- .odrl_relu_architecture(hidden)
  bounded_output <- is.character(loss) && length(loss) == 1L &&
    identical(loss, "hinge")
  score_scale <- mean(abs(score))
  if (!is.finite(score_scale)) {
    .odrl_abort("Neural network scores must be finite.")
  }
  if (score_scale <= 1e-14) {
    return(list(
      backend = "constant", backend_name = "constant-zero-score",
      constant_score = 0, architecture = architecture,
      hidden = if (length(architecture) == 1L) architecture else NA_integer_,
      activation = activation, leaky_slope = leaky_slope, decay = decay,
      loss = loss, loss_name = .odrl_relu_loss_name(loss),
      loss_builtin = is.character(loss), bounded_output = bounded_output,
      objective = 0, convergence = 0L,
      message = "Training scores were all numerically zero.",
      restart = NA_integer_, attempts = 0L, all_objectives = 0
    ))
  }
  backend <- .odrl_relu_backend_spec(backend)
  if (!backend$native) {
    model <- backend$fit(
      x = x, score = score, architecture = architecture,
      activation = activation, decay = decay, loss = loss,
      restarts = restarts, maxit = maxit, seed = seed,
      leaky_slope = leaky_slope
    )
    predicted <- backend$predict(model, x)
    if (!is.numeric(predicted) || length(predicted) != nrow(x) ||
        any(!is.finite(predicted))) {
      .odrl_abort(
        "The custom neural backend did not return one finite numeric score ",
        "per training row."
      )
    }
    model_field <- function(name, default) {
      if (is.list(model) && !is.null(model[[name]])) model[[name]] else default
    }
    model_convergence <- model_field("convergence", 0L)
    model_objective <- model_field("objective", NA_real_)
    model_restart <- model_field("restart", NA_integer_)
    model_attempts <- model_field("attempts", NA_integer_)
    model_objectives <- model_field("all_objectives", NA_real_)
    return(list(
      backend = "custom", backend_name = backend$name, model = model,
      backend_predict = backend$predict, architecture = architecture,
      hidden = if (length(architecture) == 1L) architecture else NA_integer_,
      activation = activation, leaky_slope = leaky_slope, decay = decay,
      loss = loss, loss_name = .odrl_relu_loss_name(loss),
      loss_builtin = is.character(loss), bounded_output = bounded_output,
      objective = model_objective, convergence = as.integer(model_convergence),
      message = model_field("message", NULL),
      restart = as.integer(model_restart), attempts = as.integer(model_attempts),
      all_objectives = model_objectives
    ))
  }

  # The score sign identifies the preferred treatment. Its normalized
  # magnitude determines how much that observation contributes to training.
  label <- ifelse(score >= 0, 1, -1)
  weight <- abs(score) / score_scale
  optimize_from <- function(initial, iteration_limit) {
    stats::optim(
      par = initial,
      fn = function(theta) .odrl_relu_objective(
        theta, x, label, weight, architecture, decay, loss,
        gradient = FALSE, activation = activation,
        leaky_slope = leaky_slope
      ),
      gr = function(theta) .odrl_relu_objective(
        theta, x, label, weight, architecture, decay, loss,
        gradient = TRUE, activation = activation,
        leaky_slope = leaky_slope
      ),
      method = "L-BFGS-B",
      control = list(maxit = iteration_limit, factr = 1e8)
    )
  }
  fits <- vector("list", restarts)
  for (restart in seq_len(restarts)) {
    initial <- .odrl_relu_initial(
      x, label, architecture, seed + restart * 1009L, loss, activation,
      weight = weight
    )
    fits[[restart]] <- optimize_from(initial, maxit)
  }

  objectives <- vapply(fits, `[[`, numeric(1), "value")
  converged_restarts <- which(
    vapply(fits, `[[`, integer(1), "convergence") == 0L
  )
  if (!length(converged_restarts)) {
    # If no start converged, continue the lowest-loss fit once.
    best_attempt <- which.min(objectives)
    retry <- optimize_from(
      fits[[best_attempt]]$par, max(3L * maxit, maxit + 100L)
    )
    fits[[length(fits) + 1L]] <- retry
    objectives <- c(objectives, retry$value)
    if (retry$convergence == 0L) converged_restarts <- length(fits)
  }
  # Prefer a converged fit, even if an unfinished fit has a lower objective.
  best_restart <- if (length(converged_restarts)) {
    converged_restarts[[which.min(objectives[converged_restarts])]]
  } else {
    which.min(objectives)
  }
  fit <- fits[[best_restart]]
  list(
    backend = "native", backend_name = "native", theta = fit$par,
    architecture = architecture,
    hidden = if (length(architecture) == 1L) architecture else NA_integer_,
    activation = activation, leaky_slope = leaky_slope,
    decay = decay, loss = loss, loss_name = .odrl_relu_loss_name(loss),
    loss_builtin = is.character(loss), bounded_output = bounded_output,
    objective = fit$value, convergence = fit$convergence,
    message = fit$message, restart = best_restart, attempts = length(fits),
    all_objectives = objectives
  )
}

.odrl_predict_relu_raw <- function(fit, newx) {
  raw <- if (identical(fit$backend %||% "native", "constant")) {
    rep(fit$constant_score %||% 0, nrow(newx))
  } else if (identical(fit$backend %||% "native", "custom")) {
    fit$backend_predict(fit$model, newx)
  } else {
    .odrl_relu_forward(
      fit$theta, newx, fit$architecture %||% fit$hidden,
      activation = fit$activation %||% "relu",
      leaky_slope = fit$leaky_slope %||% 0.01
    )$raw
  }
  if (!is.numeric(raw) || length(raw) != nrow(newx) ||
      any(!is.finite(raw))) {
    .odrl_abort("The neural network backend returned invalid predictions.")
  }
  if (isTRUE(fit$bounded_output)) {
    .odrl_hardtanh(raw)
  } else {
    raw
  }
}

.odrl_select_relu <- function(grid, selection) {
  converged <- which(
    grid$convergence_failures == 0L & is.finite(grid$mean_criterion)
  )
  if (!length(converged)) {
    .odrl_abort(
      "Every neural network candidate had at least one nonconverged tuning ",
      "fit. Increase `relu_maxit`, add restarts, or simplify the tuning grid."
    )
  }
  best <- converged[[which.max(grid$mean_criterion[converged])]]
  if (selection == "best") return(best)
  # The one-standard-error rule accepts candidates close to the best value,
  # then prefers fewer layers, fewer units, and stronger regularization.
  cutoff <- grid$mean_criterion[[best]] - grid$se_criterion[[best]]
  eligible <- converged[grid$mean_criterion[converged] >= cutoff]
  ordered <- eligible[order(
    grid$architecture_depth[eligible], grid$architecture_size[eligible],
    -grid$decay[eligible], grid$activation[eligible],
    -grid$mean_criterion[eligible]
  )]
  ordered[[1L]]
}

.odrl_relu_architectures_from_control <- function(control) {
  architectures <- control$relu_architectures
  if (is.null(architectures)) {
    architectures <- lapply(control$relu_hidden_units, function(width) {
      if (width == 0L) integer() else as.integer(width)
    })
  }
  if (!is.list(architectures) || !length(architectures)) {
    .odrl_abort("`relu_architectures` must be a nonempty list of layer vectors.")
  }
  lapply(architectures, .odrl_relu_architecture)
}

#' Fit a tuned neural network ODRL rule
#' @noRd
.odrl_fit_relu <- function(x, score, control, loss) {
  folds <- min(control$relu_folds, nrow(x))
  fold_id <- .odrl_score_folds(score, folds, control$seed + 400L)
  architectures <- .odrl_relu_architectures_from_control(control)
  activations <- control$relu_activation %||% "relu"
  leaky_slope <- control$relu_leaky_slope %||% 0.01
  backend <- .odrl_relu_backend_spec(
    control$relu_backend %||% "native",
    control$relu_backend_options %||% list()
  )
  grid_index <- expand.grid(
    architecture_index = seq_along(architectures),
    activation = activations,
    decay = rev(control$relu_decay),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- data.frame(
    architecture = vapply(
      architectures[grid_index$architecture_index],
      .odrl_relu_architecture_label, character(1)
    ),
    hidden = vapply(
      architectures[grid_index$architecture_index],
      function(value) if (length(value) == 1L) value else NA_integer_,
      integer(1)
    ),
    architecture_depth = lengths(
      architectures[grid_index$architecture_index]
    ),
    architecture_size = vapply(
      architectures[grid_index$architecture_index], sum, integer(1)
    ),
    activation = grid_index$activation,
    decay = grid_index$decay,
    stringsAsFactors = FALSE
  )
  grid$architecture_spec <- I(
    architectures[grid_index$architecture_index]
  )
  criterion <- matrix(NA_real_, nrow(grid), folds)
  convergence <- matrix(NA_integer_, nrow(grid), folds)
  started <- proc.time()[["elapsed"]]

  # Tune by the held-out double residual criterion. Each fold learns its own
  # input scaling and applies that transformation to the holdout rows.
  for (candidate_index in seq_len(nrow(grid))) {
    for (fold in seq_len(folds)) {
      train <- fold_id != fold
      holdout <- !train
      transform <- .odrl_minmax_fit(x[train, , drop = FALSE])
      train_x <- .odrl_minmax_apply(x[train, , drop = FALSE], transform)
      holdout_x <- .odrl_minmax_apply(x[holdout, , drop = FALSE], transform)
      fit <- .odrl_fit_relu_candidate(
        train_x, score[train], grid$architecture_spec[[candidate_index]],
        grid$decay[[candidate_index]], loss,
        control$relu_restarts, control$relu_maxit,
        control$seed + 20000L * candidate_index + fold,
        activation = grid$activation[[candidate_index]],
        leaky_slope = leaky_slope,
        backend = backend
      )
      raw <- .odrl_predict_relu_raw(fit, holdout_x)
      action <- ifelse(raw >= 0, 1, -1)
      criterion[candidate_index, fold] <- .odrl_empirical_criterion(
        score[holdout], action
      )
      convergence[candidate_index, fold] <- fit$convergence
    }
  }
  grid$mean_criterion <- rowMeans(criterion)
  grid$se_criterion <- apply(criterion, 1L, stats::sd) / sqrt(folds)
  grid$convergence_failures <- rowSums(convergence != 0L, na.rm = TRUE)
  selected <- .odrl_select_relu(grid, control$relu_selection)
  # Refit the chosen architecture and decay on all observations. The saved
  # transformation is also used for future predictions.
  transform <- .odrl_minmax_fit(x)
  transformed_x <- .odrl_minmax_apply(x, transform)
  final <- .odrl_fit_relu_candidate(
    transformed_x, score, grid$architecture_spec[[selected]],
    grid$decay[[selected]], loss, control$relu_refit_restarts,
    control$relu_maxit, control$seed + 950000L,
    activation = grid$activation[[selected]], leaky_slope = leaky_slope,
    backend = backend
  )
  if (final$convergence != 0L) {
    .odrl_abort(
      "The selected neural network candidate did not converge on the full ",
      "sample. Increase `relu_maxit` or `relu_refit_restarts`."
    )
  }
  raw <- .odrl_predict_relu_raw(final, transformed_x)
  action <- ifelse(raw >= 0, 1, -1)
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    engine = "relu",
    backend = final$backend_name,
    loss = .odrl_relu_loss_name(loss),
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
      globally_bounded = isTRUE(final$bounded_output),
      optimization_score_scale = "mean-absolute-score normalized",
      convergence = final$convergence,
      selected_restart = final$restart,
      architecture = final$architecture,
      activation = final$activation,
      leaky_slope = final$leaky_slope,
      backend = final$backend_name,
      preset = control$relu_preset %||% NA_character_,
      loss = final$loss_name
    )
  ), class = "odrl_policy_relu")
}

.odrl_predict_relu <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
  transformed_x <- .odrl_minmax_apply(newx, object$transform)
  raw <- .odrl_predict_relu_raw(object$fit, transformed_x)
  if (type == "score") raw else ifelse(raw >= 0, 1, -1)
}
