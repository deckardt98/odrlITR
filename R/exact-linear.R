.odrl_linear_design <- function(x, transform = NULL) {
  if (is.null(transform)) transform <- .odrl_standardize_fit(x)
  standardized_x <- .odrl_standardize_apply(x, transform)
  list(design = cbind(`(Intercept)` = 1, standardized_x), transform = transform)
}

.odrl_linear_solver_control <- function(control, feasibility_tolerance) {
  extra <- control$linear_solver_options %||% list()
  if (!is.list(extra) || (length(extra) &&
      (is.null(names(extra)) || any(!nzchar(names(extra)))))) {
    .odrl_abort("`linear_solver_options` must be a named list.")
  }
  reserved <- c(
    "threads", "time_limit", "log_to_console", "mip_rel_gap",
    "mip_abs_gap", "mip_feasibility_tolerance", "mip_max_nodes",
    "objective_target", "output_flag"
  )
  duplicated <- intersect(names(extra), reserved)
  if (length(duplicated)) {
    .odrl_abort(
      "`linear_solver_options` cannot override dedicated control field(s): ",
      paste(duplicated, collapse = ", "), "."
    )
  }
  arguments <- list(
    threads = control$linear_threads %||% 1L,
    time_limit = control$linear_time_limit,
    log_to_console = control$linear_log_to_console %||% FALSE,
    mip_rel_gap = control$linear_relative_gap,
    mip_feasibility_tolerance = feasibility_tolerance,
    output_flag = control$linear_log_to_console %||% FALSE
  )
  absolute_gap <- control$linear_absolute_gap
  node_limit <- control$linear_node_limit
  objective_target <- control$linear_objective_target
  if (!is.null(absolute_gap)) arguments$mip_abs_gap <- absolute_gap
  if (!is.null(node_limit)) arguments$mip_max_nodes <- node_limit
  if (!is.null(objective_target)) {
    arguments$objective_target <- objective_target
  }
  do.call(highs::highs_control, c(arguments, extra))
}

#' Fit an exact affine ODRL rule by mixed-integer optimization
#'
#' The binary decision variables encode the sign of an affine score. A finite
#' coefficient box supplies valid big-M constants, one for each row. A feasible
#' incumbent can be returned when the time limit is reached. Its solver status
#' and remaining optimization gap are recorded in `diagnostics`.
#'
#' @noRd
.odrl_fit_linear <- function(x, score, control) {
  .odrl_require("highs", "to fit an exact affine rule")
  encoded <- .odrl_linear_design(x)
  design <- encoded$design
  n <- nrow(design)
  coefficient_count <- ncol(design)
  bound <- control$linear_coefficient_bound
  margin <- control$linear_margin
  feasibility_tolerance <- min(1e-8, margin / 100)
  row_bound <- bound * rowSums(abs(design))

  # The solver vector contains coefficients followed by one binary treatment
  # indicator per row. Indicator 1 requires an affine score >= 0; indicator 0
  # requires a score <= -margin. The coefficient bounds give a valid row_bound
  # for each of these two inequalities.
  lower_constraint <- cbind(design, -diag(row_bound, nrow = n))
  upper_constraint <- cbind(design, -diag(row_bound + margin, nrow = n))
  constraint <- rbind(lower_constraint, upper_constraint)
  constraint_lower <- c(-row_bound, rep(-Inf, n))
  constraint_upper <- c(rep(Inf, n), rep(-margin, n))

  # Maximizing mean(score * indicator) gives the same rule as maximizing
  # mean(score * action), since action = 2 * indicator - 1.
  objective <- c(rep(0, coefficient_count), score / n)
  variable_lower <- c(rep(-bound, coefficient_count), rep(0, n))
  variable_upper <- c(rep(bound, coefficient_count), rep(1, n))
  variable_types <- c(rep("C", coefficient_count), rep("I", n))

  started <- proc.time()[["elapsed"]]
  solution <- highs::highs_solve(
    L = objective,
    lower = variable_lower,
    upper = variable_upper,
    A = constraint,
    lhs = constraint_lower,
    rhs = constraint_upper,
    types = variable_types,
    maximum = TRUE,
    control = .odrl_linear_solver_control(control, feasibility_tolerance)
  )
  elapsed <- proc.time()[["elapsed"]] - started
  primal <- solution$primal_solution
  primal_valid <- length(primal) == coefficient_count + n &&
    all(is.finite(primal)) &&
    isTRUE(solution$solver_msg$value_valid %||% FALSE)
  if (!primal_valid) {
    .odrl_abort(
      "The mixed-integer solver returned no finite feasible affine rule (",
      solution$status_message, ")."
    )
  }
  # Recover the rule and check that its actions match the solver's indicators.
  beta <- as.numeric(primal[seq_len(coefficient_count)])
  names(beta) <- colnames(design)
  binary_primal <- primal[coefficient_count + seq_len(n)]
  integrality_error <- max(pmin(abs(binary_primal), abs(binary_primal - 1)))
  integrality_tolerance <- max(1e-7, 10 * feasibility_tolerance)
  if (!is.finite(integrality_error) ||
      integrality_error > integrality_tolerance) {
    .odrl_abort(
      "The solver incumbent failed the binary integrality audit ",
      "(maximum distance to {0,1} = ", format(integrality_error), ")."
    )
  }
  encoded_action <- as.integer(binary_primal >= 0.5)
  affine <- drop(design %*% beta)
  action_binary <- as.integer(affine >= -feasibility_tolerance)
  constraint_value <- drop(constraint %*% primal)
  lower_violation <- max(c(constraint_lower - constraint_value, 0), na.rm = TRUE)
  upper_violation <- max(c(constraint_value - constraint_upper, 0), na.rm = TRUE)
  max_violation <- max(lower_violation, upper_violation)
  action_mismatch <- sum(action_binary != encoded_action)
  objective_audit <- mean(score * encoded_action)
  objective_difference <- abs(objective_audit - solution$objective_value)
  audit_tolerance <- max(1e-7, 10 * feasibility_tolerance)
  if (max_violation > audit_tolerance || action_mismatch > 0L ||
      objective_difference > audit_tolerance) {
    .odrl_abort(
      "The solver incumbent failed the affine rule integrity audit ",
      "(maximum constraint violation = ", format(max_violation),
      ", action mismatches = ", action_mismatch,
      ", objective difference = ", format(objective_difference), ")."
    )
  }
  action <- ifelse(action_binary == 1L, 1, -1)

  # A usable rule may still have an optimization gap. Record both the achieved
  # gap and whether it meets the user's requested stopping tolerances.
  mip_gap <- as.numeric(solution$info$mip_gap %||% NA_real_)
  solver_reported_optimal <- identical(
    tolower(solution$status_message), "optimal"
  )
  zero_gap_tolerance <- 1e-8
  certified_global_optimum <- solver_reported_optimal && is.finite(mip_gap) &&
    mip_gap <= zero_gap_tolerance
  mip_dual_bound <- as.numeric(
    solution$info$mip_dual_bound %||% NA_real_
  )
  absolute_gap <- if (is.finite(solution$objective_value) &&
      is.finite(mip_dual_bound)) {
    abs(solution$objective_value - mip_dual_bound)
  } else {
    NA_real_
  }
  requested_absolute_gap <- control$linear_absolute_gap
  relative_gap_met <- is.finite(mip_gap) &&
    mip_gap <= control$linear_relative_gap + 1e-12
  absolute_gap_met <- is.null(requested_absolute_gap) ||
    (is.finite(absolute_gap) &&
       absolute_gap <= requested_absolute_gap + 1e-12)
  requested_gap_met <- relative_gap_met && absolute_gap_met
  if (control$linear_require_gap && !requested_gap_met) {
    .odrl_abort(
      "A feasible affine rule was found, but the requested MIP gap was not ",
      "certified before the time limit. Set `linear_require_gap = FALSE` to ",
      "accept a documented incumbent."
    )
  }
  structure(list(
    engine = "highs",
    coefficients = beta,
    transform = encoded$transform,
    training_action = action,
    optimization_criterion = .odrl_empirical_criterion(score, action),
    diagnostics = list(
      solver_status = solution$status_message,
      solver_status_code = solution$status,
      solver_reported_optimal = solver_reported_optimal,
      certified_global_optimum = certified_global_optimum,
      proved_optimal = certified_global_optimum,
      zero_gap_tolerance = zero_gap_tolerance,
      requested_relative_gap = control$linear_relative_gap,
      requested_absolute_gap = requested_absolute_gap,
      mip_gap = mip_gap,
      absolute_gap = absolute_gap,
      relative_gap_met = relative_gap_met,
      absolute_gap_met = absolute_gap_met,
      requested_gap_met = requested_gap_met,
      gap_met = requested_gap_met,
      objective_value = solution$objective_value,
      mip_dual_bound = mip_dual_bound,
      mip_node_count = solution$info$mip_node_count %||% NA_real_,
      time_limit = control$linear_time_limit,
      node_limit = control$linear_node_limit,
      objective_target = control$linear_objective_target,
      threads = control$linear_threads %||% 1L,
      log_to_console = control$linear_log_to_console %||% FALSE,
      solver_options = control$linear_solver_options %||% list(),
      coefficient_bound = bound,
      margin = margin,
      feasibility_tolerance = feasibility_tolerance,
      integrality_tolerance = integrality_tolerance,
      max_integrality_error = integrality_error,
      decision_tolerance = feasibility_tolerance,
      max_constraint_violation = max_violation,
      action_mismatches = action_mismatch,
      objective_difference = objective_difference,
      optimization_score_scale = "mean-absolute-score normalized",
      elapsed = elapsed
    )
  ), class = "odrl_policy_linear")
}

.odrl_predict_linear <- function(object, newx, type = c("action", "score")) {
  type <- match.arg(type)
  design <- .odrl_linear_design(newx, object$transform)$design
  raw <- drop(design %*% object$coefficients)
  tolerance <- object$diagnostics$decision_tolerance %||% 0
  raw[abs(raw) <= tolerance] <- 0
  if (type == "score") raw else ifelse(raw >= 0, 1, -1)
}
