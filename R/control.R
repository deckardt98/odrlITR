#' Control computational settings for ODRL
#'
#' @param tree_depth Maximum depth of a direct policy tree.
#' @param tree_split_step Candidate split-step passed to
#'   [policytree::policy_tree()].
#' @param tree_min_node_size Minimum terminal-node size.
#' @param linear_coefficient_bound Absolute MILP coefficient bound.
#' @param linear_margin Numerical one-sided margin used to encode a strict
#'   negative affine score. It cannot exceed `linear_coefficient_bound`.
#' @param linear_time_limit Solver time limit in seconds.
#' @param linear_relative_gap Requested relative mixed-integer gap.
#' @param linear_require_gap If `TRUE`, reject a feasible incumbent unless the
#'   solver proves the requested gap.
#' @param svm_kernel Either `"rbf"` or `"linear"`.
#' @param svm_penalty Positive penalty grid for logistic kernel fits.
#' @param svm_rbf_multiplier Positive multipliers of the median squared
#'   pairwise distance.
#' @param svm_radius RKHS-radius grid for bounded-hinge fits. Values must not
#'   exceed one.
#' @param svm_folds Number of ODRL-criterion cross-validation folds.
#' @param svm_maxit Maximum optimizer iterations for the initial logistic
#'   kernel fit. A nonconverged fit is retried once from its incumbent with
#'   `max(3 * svm_maxit, svm_maxit + 100)` iterations.
#' @param relu_hidden_units Integer vector of one-hidden-layer widths. Include
#'   zero to add an affine score candidate.
#' @param relu_decay Nonnegative weight-decay grid.
#' @param relu_folds Number of ODRL-criterion cross-validation folds.
#' @param relu_restarts Random starts per cross-validation fit.
#' @param relu_refit_restarts Random starts for the final ReLU refit.
#' @param relu_maxit Maximum optimizer iterations per initial ReLU start. If
#'   no start converges, the best incumbent is retried once with
#'   `max(3 * relu_maxit, relu_maxit + 100)` iterations.
#' @param relu_selection Select the best criterion candidate or use a one-standard-
#'   error rule favoring smaller networks and stronger decay.
#' @param score_tolerance Scores below this absolute value are treated as zero.
#' @param seed Reproducible seed used by tuning and optimization.
#'
#' @return An object of class `odrl_control`.
#' @export
odrl_control <- function(
    tree_depth = 2L,
    tree_split_step = 1L,
    tree_min_node_size = 5L,
    linear_coefficient_bound = 10,
    linear_margin = 1e-4,
    linear_time_limit = 60,
    linear_relative_gap = 0.01,
    linear_require_gap = FALSE,
    svm_kernel = c("rbf", "linear"),
    svm_penalty = c(0.01, 0.1, 1),
    svm_rbf_multiplier = c(0.5, 1, 2),
    svm_radius = 1,
    svm_folds = 3L,
    svm_maxit = 300L,
    relu_hidden_units = c(0L, 8L, 16L),
    relu_decay = c(0.001, 0.01, 0.1),
    relu_folds = 3L,
    relu_restarts = 1L,
    relu_refit_restarts = 3L,
    relu_maxit = 250L,
    relu_selection = c("one_se", "best"),
    score_tolerance = 1e-10,
    seed = 1L) {
  svm_kernel <- match.arg(svm_kernel)
  relu_selection <- match.arg(relu_selection)
  .odrl_check_scalar(tree_depth, "tree_depth", 0, 5, integer = TRUE)
  .odrl_check_scalar(tree_split_step, "tree_split_step", 1, Inf,
                     integer = TRUE)
  .odrl_check_scalar(tree_min_node_size, "tree_min_node_size", 1, Inf,
                     integer = TRUE)
  .odrl_check_scalar(linear_coefficient_bound,
                     "linear_coefficient_bound", 0, Inf, open_lower = TRUE)
  .odrl_check_scalar(linear_margin, "linear_margin", 0, Inf,
                     open_lower = TRUE)
  if (linear_margin > linear_coefficient_bound) {
    .odrl_abort(
      "`linear_margin` cannot exceed `linear_coefficient_bound`; otherwise ",
      "even a constant negative rule may be excluded from the MILP class."
    )
  }
  .odrl_check_scalar(linear_time_limit, "linear_time_limit", 0, Inf,
                     open_lower = TRUE)
  .odrl_check_scalar(linear_relative_gap, "linear_relative_gap", 0, 1)
  if (!is.logical(linear_require_gap) || length(linear_require_gap) != 1L ||
      is.na(linear_require_gap)) {
    .odrl_abort("`linear_require_gap` must be TRUE or FALSE.")
  }
  if (!length(svm_penalty) || any(!is.finite(svm_penalty)) ||
      any(svm_penalty <= 0)) .odrl_abort("`svm_penalty` must be positive.")
  if (!length(svm_rbf_multiplier) ||
      any(!is.finite(svm_rbf_multiplier)) ||
      any(svm_rbf_multiplier <= 0)) {
    .odrl_abort("`svm_rbf_multiplier` must be positive.")
  }
  if (!length(svm_radius) || any(!is.finite(svm_radius)) ||
      any(svm_radius <= 0 | svm_radius > 1)) {
    .odrl_abort("`svm_radius` must lie in (0,1].")
  }
  .odrl_check_scalar(svm_folds, "svm_folds", 2, Inf, integer = TRUE)
  .odrl_check_scalar(svm_maxit, "svm_maxit", 1, Inf, integer = TRUE)
  if (!length(relu_hidden_units) || any(!is.finite(relu_hidden_units)) ||
      any(relu_hidden_units < 0) ||
      any(relu_hidden_units != as.integer(relu_hidden_units))) {
    .odrl_abort("`relu_hidden_units` must contain nonnegative integers.")
  }
  if (!length(relu_decay) || any(!is.finite(relu_decay)) ||
      any(relu_decay < 0)) .odrl_abort("`relu_decay` must be nonnegative.")
  .odrl_check_scalar(relu_folds, "relu_folds", 2, Inf, integer = TRUE)
  .odrl_check_scalar(relu_restarts, "relu_restarts", 1, Inf, integer = TRUE)
  .odrl_check_scalar(relu_refit_restarts, "relu_refit_restarts", 1, Inf,
                     integer = TRUE)
  .odrl_check_scalar(relu_maxit, "relu_maxit", 1, Inf, integer = TRUE)
  .odrl_check_scalar(score_tolerance, "score_tolerance", 0, Inf,
                     open_lower = TRUE)
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  structure(list(
    tree_depth = as.integer(tree_depth),
    tree_split_step = as.integer(tree_split_step),
    tree_min_node_size = as.integer(tree_min_node_size),
    linear_coefficient_bound = linear_coefficient_bound,
    linear_margin = linear_margin,
    linear_time_limit = linear_time_limit,
    linear_relative_gap = linear_relative_gap,
    linear_require_gap = linear_require_gap,
    svm_kernel = svm_kernel,
    svm_penalty = sort(unique(svm_penalty)),
    svm_rbf_multiplier = sort(unique(svm_rbf_multiplier)),
    svm_radius = sort(unique(svm_radius)),
    svm_folds = as.integer(svm_folds),
    svm_maxit = as.integer(svm_maxit),
    relu_hidden_units = sort(unique(as.integer(relu_hidden_units))),
    relu_decay = sort(unique(relu_decay)),
    relu_folds = as.integer(relu_folds),
    relu_restarts = as.integer(relu_restarts),
    relu_refit_restarts = as.integer(relu_refit_restarts),
    relu_maxit = as.integer(relu_maxit),
    relu_selection = relu_selection,
    score_tolerance = score_tolerance,
    seed = as.integer(seed)
  ), class = "odrl_control")
}

.odrl_validate_control <- function(control) {
  if (!inherits(control, "odrl_control")) {
    .odrl_abort("`control` must be created by `odrl_control()`.")
  }
  control
}
