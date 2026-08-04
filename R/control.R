#' Control computational settings for ODRL
#'
#' @param tree_depth Maximum depth of a direct policy tree.
#' @param tree_split_step Candidate split-step passed to
#'   [policytree::policy_tree()].
#' @param tree_min_node_size Minimum terminal-node size.
#' @param tree_backend Tree engine. Use `"policytree"` or a list containing
#'   `name`, `fit`, and `predict` callbacks. The fit callback receives `x`,
#'   `score`, `rewards`, `depth`, `min_node_size`, `split_step`, and `options`;
#'   the predict callback must return actions in `{-1,+1}`.
#' @param tree_options Named list of additional options for the tree engine.
#'   For `policytree`, dedicated arguments such as `depth` cannot be
#'   overridden here.
#' @param linear_coefficient_bound Absolute MILP coefficient bound.
#' @param linear_margin Numerical one-sided margin used to encode a strict
#'   negative affine score. It cannot exceed `linear_coefficient_bound`.
#' @param linear_time_limit Solver time limit in seconds.
#' @param linear_relative_gap Requested relative mixed-integer gap.
#' @param linear_require_gap If `TRUE`, reject a feasible incumbent unless the
#'   solver proves the requested gap.
#' @param linear_absolute_gap Optional requested absolute mixed-integer gap.
#' @param linear_node_limit Optional branch-and-bound node limit.
#' @param linear_objective_target Optional solver objective target.
#' @param linear_threads Positive number of HiGHS threads.
#' @param linear_log_to_console Whether HiGHS prints its solver log.
#' @param linear_solver_options Named list of additional HiGHS options. Options
#'   represented by dedicated arguments cannot be overridden here.
#' @param svm_kernel Kernel specification: `"rbf"`/`"gaussian"`, `"linear"`,
#'   `"polynomial"`/`"poly"`, a function `function(x, y)`, or a list with
#'   `name`, `fun`, and optional `args`. A custom function must return a finite
#'   kernel matrix; the caller is responsible for supplying a symmetric
#'   positive-semidefinite kernel.
#' @param svm_penalty Positive regularization grid for kernel-surrogate fits.
#' @param svm_rbf_multiplier Positive multipliers of the median squared
#'   pairwise distance.
#' @param svm_radius RKHS-radius grid for bounded-hinge fits. Values must not
#'   exceed one.
#' @param svm_folds Number of ODRL-criterion cross-validation folds.
#' @param svm_maxit Maximum optimizer iterations for the initial logistic
#'   kernel fit. A nonconverged fit is retried once from its incumbent with
#'   `max(3 * svm_maxit, svm_maxit + 100)` iterations.
#' @param svm_polynomial_degree Positive integer degree grid.
#' @param svm_polynomial_scale Positive polynomial-kernel scale.
#' @param svm_polynomial_offset Nonnegative polynomial-kernel offset.
#' @param svm_hinge_mode `"bounded"` preserves the globally bounded RBF-hinge
#'   construction; `"regularized"` enables ordinary hinge fitting with any
#'   supported kernel.
#' @param svm_kernel_function Function used with `svm_kernel = "custom"`.
#' @param svm_kernel_args Named arguments passed to a custom kernel function.
#' @param relu_hidden_units Integer vector of one-hidden-layer widths. Include
#'   zero to add an affine score candidate. This backward-compatible shortcut
#'   is ignored when `relu_architectures` is supplied.
#' @param relu_decay Nonnegative weight-decay grid.
#' @param relu_folds Number of ODRL-criterion cross-validation folds.
#' @param relu_restarts Random starts per cross-validation fit.
#' @param relu_refit_restarts Random starts for the final ReLU refit.
#' @param relu_maxit Maximum optimizer iterations per initial ReLU start. If
#'   no start converges, the best incumbent is retried once with
#'   `max(3 * relu_maxit, relu_maxit + 100)` iterations.
#' @param relu_selection Select the best criterion candidate or use a one-standard-
#'   error rule favoring smaller networks and stronger decay.
#' @param relu_architectures Optional list of hidden-layer width vectors. Use
#'   `integer()` for an affine score, `8L` for one hidden layer, and
#'   `c(16L, 8L)` for two hidden layers.
#' @param relu_activation Activation grid: `"relu"`, `"leaky_relu"`, `"tanh"`,
#'   `"sigmoid"`, or `"linear"`.
#' @param relu_leaky_slope Negative-side slope for leaky ReLU.
#' @param relu_backend `"native"` or a list containing `name`, `fit`, and
#'   `predict` callbacks for an external neural-network engine. The fit
#'   callback receives `x`, `score`, `architecture`, `activation`, `decay`,
#'   `loss`, `restarts`, `maxit`, `seed`, and `leaky_slope`. The predict
#'   callback is called as `predict(model, newx)` and must return one finite
#'   numerical score per row.
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
    seed = 1L,
    tree_backend = "policytree",
    tree_options = list(),
    linear_absolute_gap = NULL,
    linear_node_limit = NULL,
    linear_objective_target = NULL,
    linear_threads = 1L,
    linear_log_to_console = FALSE,
    linear_solver_options = list(),
    svm_polynomial_degree = 2L,
    svm_polynomial_scale = 1,
    svm_polynomial_offset = 1,
    svm_hinge_mode = c("bounded", "regularized"),
    svm_kernel_function = NULL,
    svm_kernel_args = list(),
    relu_architectures = NULL,
    relu_activation = "relu",
    relu_leaky_slope = 0.01,
    relu_backend = "native") {
  if (is.character(svm_kernel)) {
    if (!length(svm_kernel) || anyNA(svm_kernel)) {
      .odrl_abort("`svm_kernel` must contain a supported kernel name.")
    }
    svm_kernel <- tolower(svm_kernel[[1L]])
    svm_kernel <- switch(svm_kernel,
      gaussian = "rbf", poly = "polynomial", svm_kernel
    )
    if (!svm_kernel %in% c("rbf", "linear", "polynomial", "custom")) {
      .odrl_abort(
        "`svm_kernel` must be `\"rbf\"`, `\"linear\"`, `\"polynomial\"`, ",
        "`\"custom\"`, a kernel function, or a kernel specification list."
      )
    }
  } else if (!is.function(svm_kernel) && !is.list(svm_kernel)) {
    .odrl_abort(
      "`svm_kernel` must be a supported name, a function, or a list."
    )
  }
  svm_hinge_mode <- match.arg(svm_hinge_mode)
  relu_selection <- match.arg(relu_selection)
  .odrl_check_scalar(tree_depth, "tree_depth", 0, Inf, integer = TRUE)
  .odrl_check_scalar(tree_split_step, "tree_split_step", 1, Inf,
                     integer = TRUE)
  .odrl_check_scalar(tree_min_node_size, "tree_min_node_size", 1, Inf,
                     integer = TRUE)
  if (!identical(tree_backend, "policytree") &&
      (!is.list(tree_backend) || !is.function(tree_backend$fit) ||
       !is.function(tree_backend$predict))) {
    .odrl_abort(
      "`tree_backend` must be `\"policytree\"` or a list containing `fit` ",
      "and `predict` functions."
    )
  }
  if (!is.list(tree_options) || (length(tree_options) &&
      (is.null(names(tree_options)) || any(!nzchar(names(tree_options)))))) {
    .odrl_abort("`tree_options` must be a named list.")
  }
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
  if (!is.null(linear_absolute_gap)) {
    .odrl_check_scalar(linear_absolute_gap, "linear_absolute_gap", 0, Inf)
  }
  if (!is.null(linear_node_limit)) {
    .odrl_check_scalar(linear_node_limit, "linear_node_limit", 1, Inf,
                       integer = TRUE)
  }
  if (!is.null(linear_objective_target)) {
    .odrl_check_scalar(linear_objective_target, "linear_objective_target")
  }
  .odrl_check_scalar(linear_threads, "linear_threads", 1, Inf,
                     integer = TRUE)
  if (!is.logical(linear_log_to_console) ||
      length(linear_log_to_console) != 1L || is.na(linear_log_to_console)) {
    .odrl_abort("`linear_log_to_console` must be TRUE or FALSE.")
  }
  if (!is.list(linear_solver_options) || (length(linear_solver_options) &&
      (is.null(names(linear_solver_options)) ||
       any(!nzchar(names(linear_solver_options)))))) {
    .odrl_abort("`linear_solver_options` must be a named list.")
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
  if (!length(svm_polynomial_degree) ||
      any(!is.finite(svm_polynomial_degree)) ||
      any(svm_polynomial_degree < 1) ||
      any(svm_polynomial_degree != as.integer(svm_polynomial_degree))) {
    .odrl_abort("`svm_polynomial_degree` must contain positive integers.")
  }
  .odrl_check_scalar(svm_polynomial_scale, "svm_polynomial_scale", 0, Inf,
                     open_lower = TRUE)
  .odrl_check_scalar(svm_polynomial_offset, "svm_polynomial_offset", 0, Inf)
  if (!is.null(svm_kernel_function) && !is.function(svm_kernel_function)) {
    .odrl_abort("`svm_kernel_function` must be NULL or a function.")
  }
  if (!is.list(svm_kernel_args)) {
    .odrl_abort("`svm_kernel_args` must be a list.")
  }
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
  if (!is.null(relu_architectures)) {
    if (!is.list(relu_architectures) || !length(relu_architectures)) {
      .odrl_abort("`relu_architectures` must be NULL or a nonempty list.")
    }
    relu_architectures <- lapply(relu_architectures, function(architecture) {
      .odrl_relu_architecture(architecture)
    })
  }
  allowed_activation <- c("relu", "leaky_relu", "tanh", "sigmoid", "linear")
  if (!is.character(relu_activation) || !length(relu_activation) ||
      anyNA(relu_activation) ||
      any(!relu_activation %in% allowed_activation)) {
    .odrl_abort(
      "`relu_activation` must use `\"relu\"`, `\"leaky_relu\"`, `\"tanh\"`, ",
      "`\"sigmoid\"`, or `\"linear\"`."
    )
  }
  .odrl_check_scalar(relu_leaky_slope, "relu_leaky_slope", 0, 1,
                     open_lower = TRUE)
  .odrl_relu_backend_spec(relu_backend)
  .odrl_check_scalar(score_tolerance, "score_tolerance", 0, Inf,
                     open_lower = TRUE)
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  structure(list(
    tree_depth = as.integer(tree_depth),
    tree_split_step = as.integer(tree_split_step),
    tree_min_node_size = as.integer(tree_min_node_size),
    tree_backend = tree_backend,
    tree_options = tree_options,
    linear_coefficient_bound = linear_coefficient_bound,
    linear_margin = linear_margin,
    linear_time_limit = linear_time_limit,
    linear_relative_gap = linear_relative_gap,
    linear_require_gap = linear_require_gap,
    linear_absolute_gap = linear_absolute_gap,
    linear_node_limit = if (is.null(linear_node_limit)) NULL else
      as.integer(linear_node_limit),
    linear_objective_target = linear_objective_target,
    linear_threads = as.integer(linear_threads),
    linear_log_to_console = linear_log_to_console,
    linear_solver_options = linear_solver_options,
    svm_kernel = svm_kernel,
    svm_penalty = sort(unique(svm_penalty)),
    svm_rbf_multiplier = sort(unique(svm_rbf_multiplier)),
    svm_radius = sort(unique(svm_radius)),
    svm_folds = as.integer(svm_folds),
    svm_maxit = as.integer(svm_maxit),
    svm_polynomial_degree = sort(unique(as.integer(svm_polynomial_degree))),
    svm_polynomial_scale = svm_polynomial_scale,
    svm_polynomial_offset = svm_polynomial_offset,
    svm_hinge_mode = svm_hinge_mode,
    svm_kernel_function = svm_kernel_function,
    svm_kernel_args = svm_kernel_args,
    relu_hidden_units = sort(unique(as.integer(relu_hidden_units))),
    relu_decay = sort(unique(relu_decay)),
    relu_folds = as.integer(relu_folds),
    relu_restarts = as.integer(relu_restarts),
    relu_refit_restarts = as.integer(relu_refit_restarts),
    relu_maxit = as.integer(relu_maxit),
    relu_selection = relu_selection,
    relu_architectures = relu_architectures,
    relu_activation = unique(relu_activation),
    relu_leaky_slope = relu_leaky_slope,
    relu_backend = relu_backend,
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
