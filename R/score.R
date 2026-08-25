.odrl_make_score <- function(a_pm1, y, nuisance, tolerance) {
  n <- length(y)
  if (length(nuisance$m) != n || length(nuisance$e) != n) {
    .odrl_abort("Nuisance predictions must match the training sample.")
  }
  raw <- (a_pm1 - nuisance$e) * (y - nuisance$m)
  if (any(!is.finite(raw))) .odrl_abort("The ODRL score is non-finite.")
  score <- raw
  score[abs(score) < tolerance] <- 0
  scale <- mean(abs(score))
  degenerate <- !is.finite(scale) || all(score == 0)
  if (degenerate) scale <- 1
  list(
    raw = raw,
    working = score,
    scaled = score / scale,
    scale = scale,
    label = ifelse(score >= 0, 1, -1),
    weight = abs(score) / scale,
    nonzero = abs(score) >= tolerance,
    degenerate = degenerate,
    tolerance = tolerance
  )
}

#' Construct the orthogonal double residual score
#'
#' @param a Binary treatment.
#' @param y Numeric outcome.
#' @param nuisance An object from [odrl_nuisance_user()] or
#'   [odrl_nuisance_sl()].
#' @param tolerance Absolute numerical tolerance.
#' @param positive Positive treatment level for factor or character `a`.
#'
#' @return A numeric vector with the unthresholded score
#'   `(A - e_hat) * (Y - m_hat)`. Audit details are stored in attributes.
#' @export
odrl_score <- function(a, y, nuisance, tolerance = 1e-10, positive = NULL) {
  if (!inherits(nuisance, "odrl_nuisance")) {
    .odrl_abort("`nuisance` must be an `odrl_nuisance` object.")
  }
  .odrl_check_scalar(tolerance, "tolerance", 0, Inf, open_lower = TRUE)
  a_info <- .odrl_encode_treatment(a, positive = positive)
  y <- .odrl_validate_outcome(y, length(a_info$pm1))
  if (!is.null(nuisance$treatment_map) &&
      !identical(as.numeric(nuisance$treatment_map$pm1),
                 as.numeric(a_info$pm1))) {
    .odrl_abort(
      "The nuisance object was fitted under a different treatment coding ",
      "or row order. Align the positive treatment level before scoring."
    )
  }
  result <- .odrl_make_score(a_info$pm1, y, nuisance, tolerance)
  answer <- result$raw
  attr(answer, "working") <- result$working
  attr(answer, "scale") <- result$scale
  attr(answer, "degenerate") <- result$degenerate
  answer
}
