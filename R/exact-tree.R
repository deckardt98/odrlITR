#' Fit an exact ODRL policy tree
#'
#' @param x Numeric design matrix.
#' @param score Signed ODRL score.
#' @param control An [odrl_control()] object.
#'
#' @return An internal fitted-policy object.
#' @noRd
.odrl_fit_tree <- function(x, score, control) {
  .odrl_require("policytree", "to fit an exact policy tree")
  # policytree resolves equal leaf rewards in favor of the first action. Put
  # +1 first to implement the package-wide sign(0) = +1 convention.
  gamma <- cbind(`+1` = score / 2, `-1` = -score / 2)
  started <- proc.time()[["elapsed"]]
  fit <- policytree::policy_tree(
    X = x,
    Gamma = gamma,
    depth = control$tree_depth,
    split.step = control$tree_split_step,
    min.node.size = control$tree_min_node_size,
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - started
  action_id <- as.integer(stats::predict(fit, x))
  if (length(action_id) != nrow(x) || anyNA(action_id) ||
      any(!action_id %in% 1:2)) {
    .odrl_abort(
      "The policy-tree engine returned invalid action identifiers; ",
      "expected one value in {1,2} for every training row."
    )
  }
  action <- ifelse(action_id == 1L, 1, -1)
  structure(list(
    engine = "policytree",
    fit = fit,
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      depth = control$tree_depth,
      split_step = control$tree_split_step,
      min_node_size = control$tree_min_node_size,
      elapsed = elapsed,
      global_candidate_search = control$tree_split_step == 1L,
      search_scope = paste0(
        "policytree depth <= ", control$tree_depth,
        ", minimum node size ", control$tree_min_node_size,
        ", split.step ", control$tree_split_step
      )
    )
  ), class = "odrl_policy_tree")
}

.odrl_predict_tree <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
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
