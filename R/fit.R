.odrl_constant_policy <- function(n, action = 1) {
  structure(list(
    engine = "constant",
    action = action,
    training_action = rep(action, n),
    optimization_criterion = 0,
    diagnostics = list(
      reason = "All double residual scores were numerically zero.",
      tie_rule = "+1"
    )
  ), class = "odrl_policy_constant")
}

.odrl_canonical_loss_name <- function(loss) {
  name <- tolower(gsub("-", "_", loss, fixed = TRUE))
  switch(name,
    logit = "logistic",
    squaredhinge = "squared_hinge",
    sq_hinge = "squared_hinge",
    name
  )
}

.odrl_prepare_loss <- function(learner, loss) {
  if (learner %in% c("tree", "linear")) {
    if (is.null(loss)) loss <- "exact"
    if (!is.character(loss) || length(loss) != 1L || is.na(loss) ||
        !identical(tolower(loss), "exact")) {
      .odrl_abort("Direct tree and linear learners require `loss = \"exact\"`.")
    }
    return(list(fit = "exact", label = "exact"))
  }
  if (is.null(loss)) loss <- "hinge"
  if (learner == "svm") {
    specification <- .odrl_resolve_surrogate_loss(loss)
    return(list(fit = specification, label = specification$name))
  }
  if (is.character(loss)) {
    if (length(loss) != 1L || is.na(loss)) {
      .odrl_abort("`loss` must be one supported surrogate-loss name.")
    }
    name <- .odrl_canonical_loss_name(loss)
    if (!name %in% c("hinge", "exponential", "logistic", "squared_hinge")) {
      .odrl_abort(
        "Unknown neural-network surrogate loss: ", loss, ". Use `\"hinge\"`, ",
        "`\"exponential\"`, `\"logistic\"`, `\"squared_hinge\"`, or a ",
        "custom loss."
      )
    }
    return(list(fit = name, label = name))
  }
  if (is.list(loss) ||
      (is.function(loss) &&
       is.function(attr(loss, "gradient", exact = TRUE)))) {
    specification <- .odrl_resolve_surrogate_loss(loss)
    custom <- function(margin) {
      list(
        loss = specification$value(margin),
        gradient = specification$gradient(margin)
      )
    }
    attr(custom, "name") <- specification$name
    return(list(fit = custom, label = specification$name))
  }
  if (is.function(loss)) {
    .odrl_relu_custom_loss(loss, c(-1, 0, 1))
    return(list(fit = loss, label = .odrl_relu_loss_name(loss)))
  }
  .odrl_abort(
    "A custom neural-network loss must be a function or a list with ",
    "`value` and `gradient` functions."
  )
}

#' Fit an orthogonal double residual treatment rule
#'
#' `odrl()` cross-fits or accepts the propensity and marginal-outcome
#' nuisances, forms double residual scores, and fits one treatment rule.
#' Available policy classes are affine rules, shallow trees, kernel or
#' finite-series scores, and neural networks.
#'
#' @param x Covariate matrix or data frame.
#' @param a Binary treatment. Numeric values may be `{0,1}` or `{-1,+1}`.
#' @param y Numeric outcome, with larger values preferred.
#' @param learner One of `"tree"`, `"linear"`, `"svm"`, or `"relu"`.
#'   The `"relu"` learner also supports other activations and backends.
#' @param loss `"exact"` for tree/linear; `"hinge"`, `"exponential"`,
#'   `"logistic"`, or `"squared_hinge"` for SVM/neural learners; or a custom
#'   differentiable signed-margin loss. SVM custom losses use a list with
#'   `value` and `gradient` functions. Neural custom losses may use that list
#'   or a function returning `list(loss, gradient)`. If `NULL`, direct
#'   learners use exact loss and surrogate learners use hinge loss.
#' @param nuisance `NULL`/`"superlearner"` for cross-fitted Super Learner,
#'   `"parametric"`/`"glm"` for cross-fitted logistic treatment and linear
#'   outcome models, an [odrl_nuisance_user()] object, a list containing `m`
#'   plus `pi` or `e`, or a function `function(x, a, y)` returning one of those
#'   objects.
#' @param nuisance_folds Number of outer nuisance cross-fitting folds.
#' @param nuisance_fold_id Optional row-aligned outer-fold identifiers.
#'   Otherwise, treatment-stratified folds are generated.
#' @param sl.library Library passed to [odrl_nuisance_sl()].
#' @param sl.library.pi Optional propensity-specific Super Learner library.
#' @param sl.library.m Optional marginal-outcome-specific Super Learner library.
#' @param sl_inner_folds Super Learner's inner cross-validation folds.
#' @param known_pi,known_e Optional known propensity on the probability or
#'   `E(A|X)` scale. This enables randomized-trial use while cross-fitting `m`.
#' @param propensity_bounds Optional two-element limits for propensity
#'   predictions. `NULL` performs no clipping.
#' @param positive For factor/character treatment, the level representing
#'   treatment `+1`. The second factor level is used by default.
#' @param control Computational settings from [odrl_control()].
#' @param sl_verbose Whether Super Learner prints progress.
#' @param sl_env Environment for custom Super Learner wrappers.
#'
#' @return An object of class `odrl_fit`, a list with the following principal
#'   components:
#'   \describe{
#'     \item{`policy`}{The fitted second-stage policy, including its selected
#'       tuning values and learner-specific diagnostics.}
#'     \item{`nuisance`}{The aligned nuisance predictions and first-stage
#'       diagnostics used to construct the score.}
#'     \item{`score`}{Raw and normalized double residual scores, their scale,
#'       and degeneracy indicators.}
#'     \item{`training`}{Fitted treatment actions and the empirical ODRL
#'       criterion evaluated with the unscaled score.}
#'     \item{`blueprint`}{The covariate encoding retained for safe prediction.}
#'     \item{`treatment_map`}{The map between internal `{-1,+1}` actions and
#'       the treatment coding supplied by the user.}
#'   }
#'   The object also records `call`, `learner`, `loss`, `control`, sample
#'   dimensions `n` and `p`, and elapsed `runtime`. Use [predict.odrl_fit()],
#'   [fitted.odrl_fit()], [coef.odrl_fit()], and [summary.odrl_fit()] rather
#'   than relying on learner-specific internals.
#' @export
#'
#' @examples
#' data <- odrl_simulate(100, boundary = "linear", seed = 4)
#' nuisance <- odrl_nuisance_user(
#'   m = data$m, pi = data$pi, source = "known simulation truth",
#'   out_of_fold = TRUE
#' )
#' fit <- odrl(
#'   data$x, data$a, data$y, learner = "svm", nuisance = nuisance,
#'   control = odrl_control(
#'     svm_penalty = 0.1, svm_rbf_multiplier = 1, svm_folds = 2
#'   )
#' )
#' predict(fit, data$x[1:5, ])
odrl <- function(
    x, a, y,
    learner = c("tree", "linear", "svm", "relu"),
    loss = NULL,
    nuisance = NULL,
    nuisance_folds = 5L,
    sl.library = c("SL.mean", "SL.glm", "SL.glmnet"),
    sl.library.pi = sl.library,
    sl.library.m = sl.library,
    sl_inner_folds = 5L,
    known_pi = NULL,
    known_e = NULL,
    propensity_bounds = NULL,
    positive = NULL,
    control = odrl_control(),
    sl_verbose = FALSE,
    sl_env = parent.frame(),
    nuisance_fold_id = NULL) {
  call <- match.call()
  learner <- match.arg(learner)
  control <- .odrl_validate_control(control)
  loss_specification <- .odrl_prepare_loss(learner, loss)
  encoded <- .odrl_encode_x_fit(x)
  x_matrix <- encoded$x
  a_info <- .odrl_encode_treatment(a, positive = positive)
  y <- .odrl_validate_outcome(y, nrow(x_matrix))
  if (length(a_info$pm1) != nrow(x_matrix)) {
    .odrl_abort("`x`, `a`, and `y` must have matching finite rows.")
  }
  started <- proc.time()[["elapsed"]]
  nuisance_fit <- .odrl_resolve_nuisance(
    nuisance = nuisance, x = x, a = a, y = y,
    nuisance_folds = nuisance_folds, sl.library = sl.library,
    sl.library.pi = sl.library.pi, sl.library.m = sl.library.m,
    sl_inner_folds = sl_inner_folds, propensity_bounds = propensity_bounds,
    known_pi = known_pi, known_e = known_e, seed = control$seed,
    positive = positive, sl_verbose = sl_verbose, sl_env = sl_env,
    nuisance_fold_id = nuisance_fold_id
  )
  if (length(nuisance_fit$m) != nrow(x_matrix)) {
    .odrl_abort("Nuisance predictions do not match the training rows.")
  }
  if (!is.null(nuisance_fit$treatment_map) &&
      !identical(as.numeric(nuisance_fit$treatment_map$pm1),
                 as.numeric(a_info$pm1))) {
    .odrl_abort(
      "The reused nuisance object was fitted under a different treatment ",
      "coding or row order. Refit it or align the positive treatment level."
    )
  }
  if (!isTRUE(nuisance_fit$out_of_fold)) {
    warning(
      "User nuisance predictions are not marked out-of-fold or known by ",
      "design; the fitted policy does not implement the paper's cross-fitting ",
      "prescription.", call. = FALSE
    )
  }
  score <- .odrl_make_score(
    a_info$pm1, y, nuisance_fit, control$score_tolerance
  )
  policy <- if (score$degenerate) {
    .odrl_constant_policy(nrow(x_matrix), action = 1)
  } else {
    optimization_score <- score$scaled
    switch(
      learner,
      tree = .odrl_fit_tree(x_matrix, optimization_score, control),
      linear = .odrl_fit_linear(x_matrix, optimization_score, control),
      svm = .odrl_fit_svm(
        x_matrix, optimization_score, control, loss_specification$fit
      ),
      relu = .odrl_fit_relu(
        x_matrix, optimization_score, control, loss_specification$fit
      )
    )
  }
  elapsed <- proc.time()[["elapsed"]] - started
  structure(list(
    call = call,
    learner = learner,
    loss = loss_specification$label,
    policy = policy,
    nuisance = nuisance_fit,
    score = score,
    blueprint = encoded$blueprint,
    treatment_map = a_info,
    control = control,
    n = nrow(x_matrix),
    p = ncol(x_matrix),
    runtime = elapsed,
    training = list(
      action = policy$training_action,
      empirical_criterion = .odrl_empirical_criterion(
        score$raw, policy$training_action
      )
    )
  ), class = "odrl_fit")
}

#' Formula interface for ODRL
#'
#' This convenience interface supports ordinary variable-only formulas. The
#' treatment column is removed before `.` is expanded, preventing accidental
#' treatment leakage into the policy covariates.
#'
#' @param formula A two-sided formula `outcome ~ covariates`.
#' @param treatment Name of the binary treatment column.
#' @param data Data frame.
#' @param ... Arguments passed to [odrl()].
#'
#' @return An `odrl_fit` object.
#' @export
odrl_formula <- function(formula, treatment, data, ...) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    .odrl_abort("`formula` must be two-sided.")
  }
  if (length(treatment) != 1L || is.na(treatment) ||
      !nzchar(trimws(as.character(treatment)))) {
    .odrl_abort("`treatment` must be one nonempty column name.")
  }
  treatment <- as.character(treatment)
  if (!is.data.frame(data) || !treatment %in% names(data)) {
    .odrl_abort("`data` must contain the named treatment column.")
  }
  a <- data[[treatment]]
  analysis_data <- data
  analysis_data[[treatment]] <- NULL
  frame <- stats::model.frame(formula, data = analysis_data,
                              na.action = stats::na.fail)
  y <- stats::model.response(frame)
  terms <- stats::delete.response(stats::terms(frame))
  x <- stats::model.matrix(terms, data = frame)
  contrasts <- attr(x, "contrasts")
  if (attr(terms, "intercept") == 1L) {
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  }
  formula_blueprint <- list(
    kind = "formula",
    terms = terms,
    xlevels = stats::.getXlevels(terms, frame),
    contrasts = contrasts,
    columns = colnames(x),
    treatment = treatment
  )
  fit <- odrl(x = x, a = a, y = y, ...)
  fit$blueprint <- formula_blueprint
  fit$formula <- formula
  fit$treatment <- treatment
  fit
}
