#' Fit an exact ODRL policy tree
#'
#' @param x Numeric design matrix.
#' @param score Signed ODRL score.
#' @param control An [odrl_control()] object.
#'
#' @return An internal fitted-policy object.
#' @noRd
.odrl_fit_tree <- function(x, score, control) {
  # policytree resolves equal leaf rewards in favor of the first action. Put
  # +1 first to implement the package-wide sign(0) = +1 convention.
  gamma <- cbind(`+1` = score / 2, `-1` = -score / 2)
  backend <- .odrl_tree_backend_spec(control$tree_backend %||% "policytree")
  options <- control$tree_options %||% list()
  started <- proc.time()[["elapsed"]]
  if (backend$native) {
    .odrl_require("policytree", "to fit an exact policy tree")
    if (!is.list(options) || (length(options) &&
        (is.null(names(options)) || any(!nzchar(names(options)))))) {
      .odrl_abort("`tree_options` must be a named list.")
    }
    reserved <- c("X", "Gamma", "depth", "split.step", "min.node.size")
    duplicated <- intersect(names(options), reserved)
    if (length(duplicated)) {
      .odrl_abort(
        "`tree_options` cannot override dedicated control field(s): ",
        paste(duplicated, collapse = ", "), "."
      )
    }
    if (is.null(options$verbose)) options$verbose <- FALSE
    fit <- do.call(policytree::policy_tree, c(list(
      X = x,
      Gamma = gamma,
      depth = control$tree_depth,
      split.step = control$tree_split_step,
      min.node.size = control$tree_min_node_size
    ), options))
    action_id <- as.integer(stats::predict(fit, x))
    if (length(action_id) != nrow(x) || anyNA(action_id) ||
        any(!action_id %in% 1:2)) {
      .odrl_abort(
        "The policy-tree engine returned invalid action identifiers; ",
        "expected one value in {1,2} for every training row."
      )
    }
    action <- ifelse(action_id == 1L, 1, -1)
    custom_predict <- NULL
  } else {
    fit <- backend$fit(
      x = x, score = score, rewards = gamma,
      depth = control$tree_depth,
      min_node_size = control$tree_min_node_size,
      split_step = control$tree_split_step,
      options = options
    )
    action <- .odrl_tree_custom_predict(backend$predict, fit, x)
    custom_predict <- backend$predict
  }
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    engine = backend$name,
    fit = fit,
    backend_predict = custom_predict,
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      depth = control$tree_depth,
      split_step = control$tree_split_step,
      min_node_size = control$tree_min_node_size,
      backend = backend$name,
      options = options,
      elapsed = elapsed,
      global_candidate_search = backend$native &&
        control$tree_split_step == 1L,
      search_scope = paste0(
        backend$name, " depth <= ", control$tree_depth,
        ", minimum node size ", control$tree_min_node_size,
        ", split.step ", control$tree_split_step
      )
    )
  ), class = "odrl_policy_tree")
}

.odrl_tree_backend_spec <- function(backend) {
  if (is.null(backend) || identical(backend, "policytree")) {
    return(list(name = "policytree", native = TRUE))
  }
  if (!is.list(backend) || !is.function(backend$fit) ||
      !is.function(backend$predict)) {
    .odrl_abort(
      "A custom tree backend must be a list containing `fit` and `predict` ",
      "functions."
    )
  }
  name <- backend$name %||% "custom"
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    .odrl_abort("A custom tree backend's `name` must be one string.")
  }
  list(name = name, native = FALSE, fit = backend$fit,
       predict = backend$predict)
}

.odrl_tree_custom_predict <- function(predict_function, fit, newx) {
  action <- predict_function(fit, newx)
  if (!is.numeric(action) || length(action) != nrow(newx) || anyNA(action) ||
      any(!action %in% c(-1, 1))) {
    .odrl_abort(
      "A custom tree backend's `predict` function must return one action in ",
      "{-1,+1} for every row."
    )
  }
  as.numeric(action)
}

.odrl_predict_tree <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
  if (!is.null(object$backend_predict)) {
    return(.odrl_tree_custom_predict(
      object$backend_predict, object$fit, newx
    ))
  }
  action_id <- as.integer(stats::predict(object$fit, newx))
  if (length(action_id) != nrow(newx) || anyNA(action_id) ||
      any(!action_id %in% 1:2)) {
    .odrl_abort(
      "The policy-tree engine returned invalid action identifiers; ",
      "expected one value in {1,2} for every prediction row."
    )
  }
  action <- ifelse(action_id == 1L, 1, -1)
  if (type == "action") action else action
}
