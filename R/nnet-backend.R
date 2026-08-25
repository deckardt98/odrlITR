.odrl_nnet_options <- function(options = list()) {
  if (is.null(options)) options <- list()
  if (!is.list(options)) {
    .odrl_abort("`relu_backend_options` for `nnet` must be a named list.")
  }
  if (length(options) &&
      (is.null(names(options)) || any(!nzchar(names(options))))) {
    .odrl_abort("`relu_backend_options` for `nnet` must be a named list.")
  }
  defaults <- list(
    skip = TRUE,
    rang = 0.5,
    MaxNWts = 10000L,
    abstol = 1e-4,
    reltol = 1e-8,
    trace = FALSE,
    probability_epsilon = sqrt(.Machine$double.eps)
  )
  unknown <- setdiff(names(options), names(defaults))
  if (length(unknown)) {
    .odrl_abort(
      "Unknown `nnet` backend option(s): ", paste(unknown, collapse = ", "),
      ". Supported options are ", paste(names(defaults), collapse = ", "),
      "."
    )
  }
  result <- utils::modifyList(defaults, options, keep.null = FALSE)
  for (name in c("skip", "trace")) {
    value <- result[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      .odrl_abort("`nnet` option `", name, "` must be TRUE or FALSE.")
    }
  }
  .odrl_check_scalar(result$rang, "nnet option `rang`", 0, Inf)
  .odrl_check_scalar(
    result$MaxNWts, "nnet option `MaxNWts`", 1, Inf, integer = TRUE
  )
  .odrl_check_scalar(
    result$abstol, "nnet option `abstol`", 0, Inf, open_lower = TRUE
  )
  .odrl_check_scalar(
    result$reltol, "nnet option `reltol`", 0, 1, open_lower = TRUE,
    open_upper = TRUE
  )
  .odrl_check_scalar(
    result$probability_epsilon, "nnet option `probability_epsilon`",
    0, 0.5, open_lower = TRUE, open_upper = TRUE
  )
  result$MaxNWts <- as.integer(result$MaxNWts)
  result
}

.odrl_nnet_required_weights <- function(p, size, skip) {
  # Input-to-hidden and hidden-to-output links include biases. With skip
  # connections, nnet adds one direct link per input to the existing output.
  as.double((p + 1L) * size + (size + 1L) + if (skip) p else 0L)
}

.odrl_validate_nnet_candidate <- function(
    x, score, architecture, activation, decay, loss, restarts, maxit, seed,
    options) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!nrow(x) || !ncol(x) || any(!is.finite(x))) {
    .odrl_abort("The `nnet` backend requires a nonempty finite design matrix.")
  }
  if (!is.numeric(score) || length(score) != nrow(x) ||
      any(!is.finite(score))) {
    .odrl_abort(
      "The `nnet` backend requires one finite numeric score per row of `x`."
    )
  }
  if (length(architecture) &&
      (!is.numeric(architecture) || any(!is.finite(architecture)) ||
       any(architecture < 0) ||
       any(architecture != as.integer(architecture)))) {
    .odrl_abort(
      "An `nnet` architecture must be zero or one positive integer width."
    )
  }
  architecture <- .odrl_relu_architecture(architecture)
  if (length(architecture) > 1L) {
    .odrl_abort(
      "The `nnet` backend supports at most one hidden layer. Use the native ",
      "neural backend for multilayer architectures."
    )
  }
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !identical(activation, "sigmoid")) {
    .odrl_abort(
      "The `nnet` backend uses sigmoid hidden units; set ",
      "`relu_activation = \"sigmoid\"`."
    )
  }
  if (!is.character(loss) || length(loss) != 1L ||
      is.na(loss) || !identical(loss, "logistic")) {
    .odrl_abort(
      "The `nnet` backend supports only `loss = \"logistic\"`. Use the ",
      "native neural backend for hinge, exponential, squared hinge, or ",
      "custom margin losses."
    )
  }
  .odrl_check_scalar(decay, "decay", 0, Inf)
  .odrl_check_scalar(restarts, "restarts", 1, Inf, integer = TRUE)
  .odrl_check_scalar(maxit, "maxit", 1, Inf, integer = TRUE)
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  if (!length(architecture) && !isTRUE(options$skip)) {
    .odrl_abort(
      "An affine `nnet` candidate (size zero) requires `skip = TRUE`."
    )
  }
  score_scale <- mean(abs(score))
  if (!is.finite(score_scale) || score_scale <= 1e-14) {
    .odrl_abort(
      "The `nnet` backend requires at least one numerically nonzero score."
    )
  }
  size <- if (length(architecture)) architecture[[1L]] else 0L
  required_weights <- .odrl_nnet_required_weights(
    ncol(x), size, options$skip
  )
  if (required_weights > options$MaxNWts) {
    .odrl_abort(
      "The requested `nnet` candidate needs ", required_weights,
      " weights, exceeding `MaxNWts = ", options$MaxNWts, "`. Increase ",
      "`relu_backend_options$MaxNWts` or use a smaller hidden layer."
    )
  }
  list(
    x = x,
    score = as.numeric(score),
    architecture = architecture,
    size = as.integer(size),
    decay = decay,
    restarts = as.integer(restarts),
    maxit = as.integer(maxit),
    seed = as.integer(seed),
    score_scale = score_scale,
    required_weights = required_weights,
    options = options
  )
}

.odrl_nnet_restart_seed <- function(seed, restart) {
  as.integer((as.double(seed) + 1009 * as.double(restart)) %%
               .Machine$integer.max)
}

.odrl_fit_nnet_backend <- function(
    x, score, architecture, activation = "sigmoid", decay = 0,
    loss = "logistic", restarts = 1L, maxit = 100L, seed = 1L,
    leaky_slope = 0.01, options = list()) {
  .odrl_require("nnet", "to use the `nnet` neural backend")
  options <- .odrl_nnet_options(options)
  specification <- .odrl_validate_nnet_candidate(
    x = x, score = score, architecture = architecture,
    activation = activation, decay = decay, loss = loss,
    restarts = restarts, maxit = maxit, seed = seed, options = options
  )
  # The score is converted to weighted binary classification. Scaling the
  # case weights to sum approximately one makes nnet's entropy term equal to
  # the package's mean normalized weighted logistic margin loss.
  label <- as.numeric(specification$score >= 0)
  case_weight <- abs(specification$score) /
    (specification$score_scale * nrow(specification$x))
  fit_once <- function(restart, maxit, Wts = NULL) {
    set.seed(.odrl_nnet_restart_seed(specification$seed, restart))
    arguments <- list(
      x = specification$x,
      y = label,
      weights = case_weight,
      size = specification$size,
      entropy = TRUE,
      linout = FALSE,
      softmax = FALSE,
      skip = options$skip,
      rang = options$rang,
      decay = specification$decay,
      maxit = maxit,
      Hess = FALSE,
      trace = options$trace,
      MaxNWts = options$MaxNWts,
      abstol = options$abstol,
      reltol = options$reltol
    )
    if (!is.null(Wts)) arguments$Wts <- Wts
    do.call(nnet::nnet, arguments)
  }
  fits <- lapply(seq_len(specification$restarts), function(restart) {
    fit_once(restart, specification$maxit)
  })
  values <- vapply(fits, `[[`, numeric(1), "value")
  convergence <- vapply(fits, `[[`, integer(1), "convergence")
  converged <- which(convergence == 0L & is.finite(values))
  if (!length(converged)) {
    finite <- which(is.finite(values))
    if (!length(finite)) {
      .odrl_abort("Every `nnet` restart returned a nonfinite objective.")
    }
    provisional <- finite[[which.min(values[finite])]]
    retry <- fit_once(
      specification$restarts + 1L,
      max(3L * specification$maxit, specification$maxit + 100L),
      Wts = fits[[provisional]]$wts
    )
    fits[[length(fits) + 1L]] <- retry
    values <- c(values, retry$value)
    convergence <- c(convergence, retry$convergence)
    if (retry$convergence == 0L && is.finite(retry$value)) {
      converged <- length(fits)
    }
  }
  eligible <- if (length(converged)) converged else which(is.finite(values))
  if (!length(eligible)) {
    .odrl_abort("Every `nnet` restart returned a nonfinite objective.")
  }
  best <- eligible[[which.min(values[eligible])]]
  structure(list(
    fit = fits[[best]],
    architecture = specification$architecture,
    activation = "sigmoid",
    decay = specification$decay,
    loss = "logistic",
    convergence = as.integer(convergence[[best]]),
    objective = values[[best]],
    restart = as.integer(best),
    attempts = length(fits),
    all_objectives = values,
    all_convergence = as.integer(convergence),
    score_scale = specification$score_scale,
    case_weight_sum = sum(case_weight),
    required_weights = specification$required_weights,
    options = options,
    leaky_slope = leaky_slope
  ), class = "odrl_nnet_backend_fit")
}

.odrl_predict_nnet_backend <- function(model, newx) {
  if (!inherits(model, "odrl_nnet_backend_fit") ||
      !inherits(model$fit, "nnet")) {
    .odrl_abort("`model` is not a fitted odrlITR `nnet` backend object.")
  }
  newx <- as.matrix(newx)
  storage.mode(newx) <- "double"
  if (any(!is.finite(newx))) {
    .odrl_abort("`newx` for the `nnet` backend must be finite.")
  }
  probability <- drop(stats::predict(model$fit, newx, type = "raw"))
  if (!is.numeric(probability) || length(probability) != nrow(newx) ||
      any(!is.finite(probability))) {
    .odrl_abort(
      "The `nnet` backend did not return one finite probability per row."
    )
  }
  epsilon <- model$options$probability_epsilon
  stats::qlogis(pmin(pmax(probability, epsilon), 1 - epsilon))
}

.odrl_nnet_backend <- function(options = list()) {
  options <- .odrl_nnet_options(options)
  list(
    name = "nnet",
    fit = function(x, score, architecture, activation, decay, loss,
                   restarts, maxit, seed, leaky_slope) {
      .odrl_fit_nnet_backend(
        x = x, score = score, architecture = architecture,
        activation = activation, decay = decay, loss = loss,
        restarts = restarts, maxit = maxit, seed = seed,
        leaky_slope = leaky_slope, options = options
      )
    },
    predict = function(model, newx) {
      .odrl_predict_nnet_backend(model, newx)
    }
  )
}
