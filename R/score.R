.odrl_make_score <- function(a_pm1, y, nuisance, tolerance) {
  n <- length(y)
  if (length(nuisance$m) != n || length(nuisance$e) != n) {
    .odrl_abort("Nuisance predictions must match the training sample.")
  }
  # Remove the predictable part of treatment and outcome, then multiply the
  # residuals. Keep this raw score for reporting the empirical criterion.
  treatment_residual <- a_pm1 - nuisance$e
  outcome_residual <- y - nuisance$m
  raw_score <- treatment_residual * outcome_residual
  if (any(!is.finite(raw_score))) {
    .odrl_abort("The ODRL score is non-finite.")
  }

  # Ignore numerical zeros during fitting, then scale by the mean absolute
  # score across all rows. Signs become labels and magnitudes become weights.
  working_score <- raw_score
  working_score[abs(working_score) < tolerance] <- 0
  scale <- mean(abs(working_score))
  degenerate <- !is.finite(scale) || all(working_score == 0)
  if (degenerate) scale <- 1

  list(
    raw = raw_score,
    working = working_score,
    scaled = working_score / scale,
    scale = scale,
    label = ifelse(working_score >= 0, 1, -1),
    weight = abs(working_score) / scale,
    nonzero = abs(working_score) >= tolerance,
    degenerate = degenerate,
    tolerance = tolerance
  )
}

#' Construct the orthogonal double residual score
#'
#' The score is `(A - e_hat) * (Y - m_hat)`, where treatment `A` is coded
#' as `-1` or `+1`, `e_hat` predicts `E(A | X)`, and `m_hat` predicts
#' `E(Y | X)`. The sign determines the treatment label used for fitting,
#' and the absolute value determines its weight.
#'
#' @param a Binary treatment.
#' @param y Numeric outcome.
#' @param nuisance An object from [odrl_nuisance_user()],
#'   [odrl_nuisance_sl()], or [odrl_nuisance_parametric()]. Predictions must
#'   follow the order of `a` and `y`.
#' @param tolerance Positive threshold below which scores are treated as zero
#'   during fitting. The returned score retains the original values.
#' @param positive Positive treatment level for factor or character `a`.
#'
#' @return A numeric vector of raw scores. Its `working` attribute contains
#'   the scores after applying `tolerance`, `scale` is their mean absolute
#'   value, and `degenerate` indicates that all working scores are zero.
#'   A degenerate score uses `scale = 1`.
#' @export
odrl_score <- function(a, y, nuisance, tolerance = 1e-10, positive = NULL) {
  if (!inherits(nuisance, "odrl_nuisance")) {
    .odrl_abort("`nuisance` must be an `odrl_nuisance` object.")
  }
  .odrl_check_scalar(tolerance, "tolerance", 0, Inf, open_lower = TRUE)
  treatment_info <- .odrl_encode_treatment(a, positive = positive)
  y <- .odrl_validate_outcome(y, length(treatment_info$pm1))
  if (!is.null(nuisance$treatment_map) &&
      !identical(as.numeric(nuisance$treatment_map$pm1),
                 as.numeric(treatment_info$pm1))) {
    .odrl_abort(
      "The nuisance object was fitted under a different treatment coding ",
      "or row order. Align the positive treatment level before scoring."
    )
  }
  result <- .odrl_make_score(treatment_info$pm1, y, nuisance, tolerance)
  answer <- result$raw
  attr(answer, "working") <- result$working
  attr(answer, "scale") <- result$scale
  attr(answer, "degenerate") <- result$degenerate
  answer
}
