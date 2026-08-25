#' Predict from an ODRL fit
#'
#' `type = "action"` returns treatment recommendations in the coding used to
#' fit the model. `type = "score"` returns the real-valued decision score for
#' affine and surrogate learners. Trees and constant policies have no unique
#' real-valued score, so their signed action is returned. For hinge SVMs with
#' bounded scores, the returned scores are clipped to `[-1,1]` with hard tanh.
#' A zero score maps to the positive treatment.
#'
#' @param object An [odrl()] fit.
#' @param newdata New covariates in the same representation used for fitting.
#' @param type Return decoded treatment actions or decision scores.
#' @param ... Unused.
#'
#' @return A vector of treatment actions or numeric decision scores.
#' @export
predict.odrl_fit <- function(object, newdata, type = c("action", "score"),
                             ...) {
  type <- match.arg(type)
  newx <- .odrl_encode_x_new(newdata, object$blueprint)
  if (inherits(object$policy, "odrl_policy_constant")) {
    raw <- if (type == "score") {
      rep(object$policy$action, nrow(newx))
    } else {
      rep(object$policy$action, nrow(newx))
    }
  } else {
    raw <- switch(
      object$learner,
      tree = .odrl_predict_tree(object$policy, newx, type = type),
      linear = .odrl_predict_linear(object$policy, newx, type = type),
      svm = .odrl_predict_svm(object$policy, newx, type = type),
      relu = .odrl_predict_relu(object$policy, newx, type = type)
    )
  }
  if (type == "score") return(as.numeric(raw))
  .odrl_decode_action(as.numeric(raw), object$treatment_map)
}

#' Extract fitted treatment recommendations
#'
#' @param object An [odrl()] fit.
#' @param ... Unused.
#'
#' @return Treatment recommendations in the original treatment coding.
#' @export
fitted.odrl_fit <- function(object, ...) {
  .odrl_decode_action(object$policy$training_action, object$treatment_map)
}

#' Extract affine rule coefficients
#'
#' By default, coefficients are returned on the original encoded covariate
#' scale. Use `standardized = TRUE` to obtain the coefficients used internally
#' after standardization.
#'
#' @param object An [odrl()] fit.
#' @param standardized Whether to return internal standardized coefficients.
#' @param ... Unused.
#'
#' @return A named numeric vector for affine fits, otherwise `NULL`.
#' @export
coef.odrl_fit <- function(object, standardized = FALSE, ...) {
  if (!inherits(object$policy, "odrl_policy_linear")) return(NULL)
  beta <- object$policy$coefficients
  if (standardized) return(beta)
  transform <- object$policy$transform
  slope <- beta[-1L] / transform$scale
  intercept <- beta[[1L]] - sum(beta[-1L] * transform$center /
                                  transform$scale)
  answer <- c(`(Intercept)` = intercept, slope)
  names(answer) <- names(beta)
  answer
}

#' Print an ODRL fit
#'
#' @param x An [odrl()] fit.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.odrl_fit <- function(x, ...) {
  cat("Orthogonal double residual learning fit\n")
  cat("  learner:  ", x$learner, "\n", sep = "")
  cat("  loss:     ", x$loss, "\n", sep = "")
  cat("  sample:   ", x$n, " rows; ", x$p, " encoded covariates\n", sep = "")
  cat("  nuisance: ", x$nuisance$source, "\n", sep = "")
  cat("  criterion:", format(x$training$empirical_criterion, digits = 5), "\n",
      sep = "")
  cat("  runtime:  ", .odrl_format_seconds(x$runtime), "\n", sep = "")
  invisible(x)
}

#' Summarize an ODRL fit
#'
#' @param object An [odrl()] fit.
#' @param ... Unused.
#' @return An object of class `summary.odrl_fit`.
#' @export
summary.odrl_fit <- function(object, ...) {
  structure(list(
    learner = object$learner,
    loss = object$loss,
    n = object$n,
    p = object$p,
    nuisance_source = object$nuisance$source,
    nuisance_out_of_fold = object$nuisance$out_of_fold,
    propensity_range = object$nuisance$diagnostics$pi_range,
    clipped_fraction = object$nuisance$diagnostics$clipped_fraction,
    score_scale = object$score$scale,
    score_zero_fraction = mean(!object$score$nonzero),
    empirical_criterion = object$training$empirical_criterion,
    selected_tuning = object$policy$selected %||% NULL,
    policy_diagnostics = object$policy$diagnostics,
    runtime = object$runtime
  ), class = "summary.odrl_fit")
}

#' @export
print.summary.odrl_fit <- function(x, ...) {
  cat("ODRL fit summary\n")
  cat("  learner/loss: ", x$learner, " / ", x$loss, "\n", sep = "")
  cat("  n / p:        ", x$n, " / ", x$p, "\n", sep = "")
  cat("  nuisance:     ", x$nuisance_source, "\n", sep = "")
  cat("  out-of-fold:   ", x$nuisance_out_of_fold, "\n", sep = "")
  cat("  propensity:    [", paste(format(x$propensity_range, digits = 4),
                                  collapse = ", "), "]\n", sep = "")
  cat("  ODRL criterion: ", format(x$empirical_criterion, digits = 5), "\n",
      sep = "")
  if (!is.null(x$selected_tuning)) {
    cat("  selected tuning:\n")
    print(x$selected_tuning, row.names = FALSE)
  }
  cat("  runtime:       ", .odrl_format_seconds(x$runtime), "\n", sep = "")
  invisible(x)
}
