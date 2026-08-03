.odrl_linear_design <- function(x, transform = NULL) {
  if (is.null(transform)) transform <- .odrl_standardize_fit(x)
  z <- .odrl_standardize_apply(x, transform)
  list(design = cbind(`(Intercept)` = 1, z), transform = transform)
}

#' Fit an exact affine ODRL rule by mixed-integer optimization
#'
#' The binary decision variables encode the sign of an affine score. A finite
#' coefficient box supplies valid row-specific big-M constants. A feasible
#' incumbent is returned when the time limit is reached; its solver status and
#' MIP gap remain available in `diagnostics` and are never relabelled optimal.
#'
#' @noRd
.odrl_fit_linear <- function(x, score, control) {
  .odrl_require("highs", "to fit an exact affine rule")
  encoded <- .odrl_linear_design(x)
  b <- encoded$design
  n <- nrow(b)
  q <- ncol(b)
  bound <- control$linear_coefficient_bound
  margin <- control$linear_margin
  feasibility_tolerance <- min(1e-8, margin / 100)
  row_bound <- bound * rowSums(abs(b))

  lower_constraint <- cbind(b, -diag(row_bound, nrow = n))
  upper_constraint <- cbind(b, -diag(row_bound + margin, nrow = n))
  constraint <- rbind(lower_constraint, upper_constraint)
  lhs <- c(-row_bound, rep(-Inf, n))
  rhs <- c(rep(Inf, n), rep(-margin, n))
  objective <- c(rep(0, q), score / n)
  lower <- c(rep(-bound, q), rep(0, n))
  upper <- c(rep(bound, q), rep(1, n))
  types <- c(rep("C", q), rep("I", n))

  started <- proc.time()[["elapsed"]]
  solution <- highs::highs_solve(
    L = objective,
    lower = lower,
    upper = upper,
    A = constraint,
    lhs = lhs,
    rhs = rhs,
    types = types,
    maximum = TRUE,
    control = highs::highs_control(
      threads = 1L,
      time_limit = control$linear_time_limit,
      log_to_console = FALSE,
      mip_rel_gap = control$linear_relative_gap,
      mip_feasibility_tolerance = feasibility_tolerance,
      output_flag = FALSE
    )
  )
  elapsed <- proc.time()[["elapsed"]] - started
  primal <- solution$primal_solution
  primal_valid <- length(primal) == q + n && all(is.finite(primal)) &&
    isTRUE(solution$solver_msg$value_valid %||% FALSE)
  if (!primal_valid) {
    .odrl_abort(
      "The mixed-integer solver returned no finite feasible affine rule (",
      solution$status_message, ")."
    )
  }
  beta <- as.numeric(primal[seq_len(q)])
  names(beta) <- colnames(b)
  binary_primal <- primal[q + seq_len(n)]
  integrality_error <- max(pmin(abs(binary_primal), abs(binary_primal - 1)))
  integrality_tolerance <- max(1e-7, 10 * feasibility_tolerance)
  if (!is.finite(integrality_error) ||
      integrality_error > integrality_tolerance) {
    .odrl_abort(
      "The solver incumbent failed the binary-integrality audit ",
      "(maximum distance to {0,1} = ", format(integrality_error), ")."
    )
  }
  encoded_action <- as.integer(binary_primal >= 0.5)
  affine <- drop(b %*% beta)
  action_binary <- as.integer(affine >= -feasibility_tolerance)
  constraint_value <- drop(constraint %*% primal)
  lower_violation <- max(c(lhs - constraint_value, 0), na.rm = TRUE)
  upper_violation <- max(c(constraint_value - rhs, 0), na.rm = TRUE)
  max_violation <- max(lower_violation, upper_violation)
  action_mismatch <- sum(action_binary != encoded_action)
  objective_audit <- mean(score * encoded_action)
  objective_difference <- abs(objective_audit - solution$objective_value)
  audit_tolerance <- max(1e-7, 10 * feasibility_tolerance)
  if (max_violation > audit_tolerance || action_mismatch > 0L ||
      objective_difference > audit_tolerance) {
    .odrl_abort(
      "The solver incumbent failed the affine-rule integrity audit ",
      "(maximum constraint violation = ", format(max_violation),
      ", action mismatches = ", action_mismatch,
      ", objective difference = ", format(objective_difference), ")."
    )
  }
  action <- ifelse(action_binary == 1L, 1, -1)
  mip_gap <- as.numeric(solution$info$mip_gap %||% NA_real_)
  solver_reported_optimal <- identical(
    tolower(solution$status_message), "optimal"
  )
  zero_gap_tolerance <- 1e-8
  certified_global_optimum <- solver_reported_optimal && is.finite(mip_gap) &&
    mip_gap <= zero_gap_tolerance
  requested_gap_met <- is.finite(mip_gap) &&
    mip_gap <= control$linear_relative_gap + 1e-12
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
      mip_gap = mip_gap,
      requested_gap_met = requested_gap_met,
      gap_met = requested_gap_met,
      objective_value = solution$objective_value,
      mip_dual_bound = solution$info$mip_dual_bound %||% NA_real_,
      mip_node_count = solution$info$mip_node_count %||% NA_real_,
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
