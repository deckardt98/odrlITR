#' Supply treatment and outcome predictions to ODRL
#'
#' Use this function when you have already obtained the treatment probability
#' and mean outcome predictions. Supply one mean outcome prediction `m` per
#' row, in the same order as the data passed to [odrl()]. For each estimated
#' quantity, predict a row using a model fitted without that row. Known
#' quantities, such as a randomized treatment probability, can be supplied
#' directly.
#'
#' Supply either the treatment probability `pi` or the expected treatment `e`.
#' These are related by `e = 2 * pi - 1` because the package uses treatment
#' values \eqn{\{-1,+1\}} internally.
#'
#' @param m Numeric predictions of \eqn{E(Y\mid X)}.
#' @param pi Numeric predictions of \eqn{P(A=+1\mid X)}.
#' @param e Numeric predictions of \eqn{E(A\mid X)}.
#' @param fold_id Optional identifier of the held-out fold for each row.
#' @param source Short description of how the predictions were obtained.
#' @param out_of_fold Set to `TRUE` when every supplied prediction was
#'   obtained from a model fitted without that row or is known by design.
#'   The default is `FALSE`, which causes [odrl()] to warn that cross-fitting
#'   has not been confirmed.
#' @param propensity_bounds Optional numeric vector of length two. Use
#'   `NULL` to retain supplied probabilities exactly.
#'
#' @return An `odrl_nuisance` object containing `m`, `pi`, and `e`, along with
#'   the prediction source, fold identifiers, and clipping diagnostics. A
#'   scalar `pi` or `e` is repeated for all rows.
#' @export
odrl_nuisance_user <- function(
    m, pi = NULL, e = NULL, fold_id = NULL,
    source = "user supplied", out_of_fold = FALSE,
    propensity_bounds = NULL) {
  if (!is.numeric(m) || is.factor(m)) {
    .odrl_abort("`m` must be a genuinely numeric prediction vector.")
  }
  m <- as.numeric(m)
  if (!length(m) || any(!is.finite(m))) {
    .odrl_abort("`m` must be a nonempty finite numeric vector.")
  }
  if (is.null(pi) == is.null(e)) {
    .odrl_abort("Supply exactly one of `pi` and `e`.")
  }

  # Convert either treatment scale to probabilities before applying bounds.
  if (!is.null(pi)) {
    if (!is.numeric(pi) || is.factor(pi)) {
      .odrl_abort("`pi` must be a genuinely numeric prediction vector.")
    }
    pi <- as.numeric(pi)
    if (length(pi) == 1L) pi <- rep(pi, length(m))
    if (length(pi) != length(m) || any(!is.finite(pi)) ||
        any(pi < 0 | pi > 1)) {
      .odrl_abort("`pi` must match `m` and lie in [0,1].")
    }
    raw_pi <- pi
  } else {
    if (!is.numeric(e) || is.factor(e)) {
      .odrl_abort("`e` must be a genuinely numeric prediction vector.")
    }
    e <- as.numeric(e)
    if (length(e) == 1L) e <- rep(e, length(m))
    if (length(e) != length(m) || any(!is.finite(e)) ||
        any(e < -1 | e > 1)) {
      .odrl_abort("`e` must match `m` and lie in [-1,1].")
    }
    raw_pi <- (e + 1) / 2
  }

  # Preserve raw probabilities so the effect of optional clipping is visible.
  if (!is.null(propensity_bounds)) {
    if (length(propensity_bounds) != 2L ||
        any(!is.finite(propensity_bounds)) ||
        propensity_bounds[[1L]] < 0 || propensity_bounds[[2L]] > 1 ||
        propensity_bounds[[1L]] >= propensity_bounds[[2L]]) {
      .odrl_abort("`propensity_bounds` must be increasing inside [0,1].")
    }
    pi <- pmin(pmax(raw_pi, propensity_bounds[[1L]]),
               propensity_bounds[[2L]])
  } else {
    pi <- raw_pi
  }
  if (!is.null(fold_id)) {
    .odrl_validate_fold_id(fold_id, length(m), normalize = FALSE)
  }
  if (!is.logical(out_of_fold) || length(out_of_fold) != 1L ||
      is.na(out_of_fold)) {
    .odrl_abort("`out_of_fold` must be TRUE or FALSE.")
  }
  if (length(source) != 1L || is.na(source) ||
      !nzchar(trimws(as.character(source)))) {
    .odrl_abort("`source` must be one nonempty description.")
  }
  structure(list(
    m = m,
    pi = pi,
    e = 2 * pi - 1,
    raw_pi = raw_pi,
    fold_id = fold_id,
    source = as.character(source),
    out_of_fold = out_of_fold,
    propensity_bounds = propensity_bounds,
    diagnostics = list(
      pi_range = range(pi),
      clipped_fraction = mean(pi != raw_pi)
    )
  ), class = "odrl_nuisance")
}

.odrl_validate_fold_id <- function(
    fold_id, n, a_pm1 = NULL, require_training_arms = FALSE,
    min_training_n = 0L, normalize = TRUE) {
  valid_type <- is.numeric(fold_id) || is.character(fold_id) ||
    is.factor(fold_id) || is.logical(fold_id)
  if (is.data.frame(fold_id) || is.list(fold_id) || !valid_type) {
    .odrl_abort("`fold_id` must be an atomic vector with one value per row.")
  }
  if (length(fold_id) != n) {
    .odrl_abort("`fold_id` must match the number of rows.")
  }
  if (anyNA(fold_id) ||
      (is.numeric(fold_id) && any(!is.finite(fold_id)))) {
    .odrl_abort("`fold_id` cannot contain missing or non-finite values.")
  }
  row_labels <- as.character(fold_id)
  if (any(!nzchar(trimws(row_labels)))) {
    .odrl_abort("`fold_id` cannot contain empty labels.")
  }
  fold_labels <- unique(row_labels)
  if (length(fold_labels) < 2L) {
    .odrl_abort("`fold_id` must identify at least two nonempty folds.")
  }

  # Normalize labels in order of first appearance, preserving that order in
  # the diagnostics even when the supplied labels are factors.
  id <- match(row_labels, fold_labels)
  if (min_training_n > 0L) {
    training_sizes <- vapply(
      seq_along(fold_labels), function(fold) sum(id != fold), integer(1)
    )
    if (any(training_sizes < min_training_n)) {
      minimum_label <- if (identical(as.integer(min_training_n), 2L)) {
        "two"
      } else {
        as.character(min_training_n)
      }
      .odrl_abort(
        "Each user-specified outer training set must contain at least ",
        minimum_label, " observations."
      )
    }
  }
  if (!is.null(a_pm1) && require_training_arms) {
    for (fold in seq_along(fold_labels)) {
      arm_counts <- table(factor(a_pm1[id != fold], levels = c(-1, 1)))
      if (min(arm_counts) < 2L) {
        .odrl_abort(
          "Each user-specified outer training set must contain at least two ",
          "observations from each treatment arm."
        )
      }
    }
  }
  list(
    id = if (normalize) as.integer(id) else fold_id,
    labels = fold_labels,
    folds = length(fold_labels),
    source = "user supplied"
  )
}

.odrl_treatment_stratified_folds <- function(
    a_pm1, folds, seed, require_training_arms = TRUE) {
  folds <- as.integer(folds)
  .odrl_check_scalar(folds, "folds", 2, length(a_pm1), integer = TRUE)
  arm_counts <- table(factor(a_pm1, levels = c(-1, 1)))
  if (require_training_arms && min(arm_counts) < folds) {
    .odrl_abort("Each treatment arm must contain at least `folds` rows.")
  }

  set.seed(seed)
  fold_id <- integer(length(a_pm1))
  fold_sizes <- integer(folds)
  for (arm in c(-1, 1)) {
    arm_rows <- sample(which(a_pm1 == arm))
    rows_per_fold <- length(arm_rows) %/% folds
    extra_rows <- length(arm_rows) %% folds
    arm_fold_sizes <- rep.int(rows_per_fold, folds)
    if (extra_rows > 0L) {
      # Spread each arm evenly, giving its extra rows to the smallest folds.
      # Random tie breaking keeps equally sized folds interchangeable.
      priority <- order(fold_sizes, stats::runif(folds))
      smallest_folds <- priority[seq_len(extra_rows)]
      arm_fold_sizes[smallest_folds] <- arm_fold_sizes[smallest_folds] + 1L
    }
    fold_labels <- rep.int(seq_len(folds), arm_fold_sizes)
    fold_id[arm_rows] <- sample(
      fold_labels, length(fold_labels), replace = FALSE
    )
    fold_sizes <- fold_sizes + arm_fold_sizes
  }
  fold_id
}

.odrl_crossfit_fold_spec <- function(a_pm1, folds, fold_id, seed,
                                     require_training_arms = TRUE) {
  if (!is.null(fold_id)) {
    return(.odrl_validate_fold_id(
      fold_id, length(a_pm1), a_pm1 = a_pm1,
      require_training_arms = require_training_arms,
      min_training_n = 2L, normalize = TRUE
    ))
  }

  # Each fold supplies predictions only for its held-out rows. Check the
  # complementary training sets before any model fitting starts.
  id <- .odrl_treatment_stratified_folds(
    a_pm1, folds, seed,
    require_training_arms = require_training_arms
  )
  fold_numbers <- seq_len(as.integer(folds))
  training_sizes <- vapply(
    fold_numbers, function(fold) sum(id != fold), integer(1)
  )
  if (any(training_sizes < 2L)) {
    .odrl_abort(
      "Each outer training set must contain at least two observations; ",
      "reduce `folds`."
    )
  }
  if (require_training_arms) {
    for (fold in fold_numbers) {
      arm_counts <- table(factor(a_pm1[id != fold], levels = c(-1, 1)))
      if (min(arm_counts) < 2L) {
        .odrl_abort(
          "Each outer training set must contain at least two observations ",
          "from each treatment arm; reduce `folds`."
        )
      }
    }
  }
  list(
    id = id,
    labels = as.character(fold_numbers),
    folds = as.integer(folds),
    source = "treatment-stratified"
  )
}

# Convert a known treatment probability or mean to one probability per row.
# Both nuisance engines use this conversion and skip propensity fitting when
# it returns a vector. NULL means that the propensity must be estimated.
.odrl_known_propensity <- function(known_pi, known_e, n) {
  if (!is.null(known_pi) && !is.null(known_e)) {
    .odrl_abort("Supply at most one of `known_pi` and `known_e`.")
  }
  if (is.null(known_pi) && is.null(known_e)) return(NULL)

  supplied <- if (!is.null(known_pi)) known_pi else known_e
  if (!is.numeric(supplied) || is.factor(supplied)) {
    .odrl_abort("Known propensity values must be genuinely numeric.")
  }
  probability <- if (!is.null(known_pi)) {
    as.numeric(known_pi)
  } else {
    (as.numeric(known_e) + 1) / 2
  }
  if (length(probability) == 1L) probability <- rep(probability, n)
  if (length(probability) != n || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    .odrl_abort(
      "Known propensity values must align with rows and lie in [0,1]."
    )
  }
  probability
}

#' Estimate treatment and outcome predictions with Super Learner
#'
#' Fits the treatment probability \eqn{P(A=+1\mid X)} and mean outcome
#' \eqn{E(Y\mid X)} using cross-fitting. Each fold is held out in turn while
#' models are fitted on the remaining rows. The resulting predictions follow
#' the original row order and can be passed directly to [odrl()]. The outcome
#' model pools both treatment arms and uses covariates alone.
#'
#' Within each training set, Super Learner uses its own cross-validation
#' folds to combine the candidate models in its library. `folds` controls
#' the outer prediction folds, and `inner_folds` controls this model
#' combination step.
#'
#' @param x Covariate matrix or data frame.
#' @param a Binary treatment, coded as `{0,1}`, `{-1,+1}`, logical, or a
#'   two-level factor.
#' @param y Numeric outcome.
#' @param folds Number of outer cross-fitting folds. Ignored when `fold_id` is
#'   supplied. Automatically generated folds are treatment-stratified, and
#'   their total sizes differ by at most one row.
#' @param fold_id Optional outer fold identifiers supplied for each row.
#'   Labels may be numeric, character, or factor and are normalized internally.
#'   Every resulting training set must contain at least two observations and
#'   must retain at least two observations from each treatment arm when the
#'   propensity is estimated.
#' @param sl.library Character vector or list understood by
#'   [SuperLearner::SuperLearner()].
#' @param sl.library.pi Optional library for the propensity. Defaults to
#'   `sl.library`.
#' @param sl.library.m Optional library for the marginal outcome. Defaults to
#'   `sl.library`.
#' @param inner_folds Super Learner cross-validation folds within each outer
#'   training sample.
#' @param seed Reproducibility seed.
#' @param positive For a factor or character treatment, the level representing
#'   treatment `+1`. By default, the second factor level is positive.
#' @param propensity_bounds Optional lower and upper limits for predicted
#'   treatment probabilities, such as `c(0.01, 0.99)`. The default `NULL`
#'   keeps the predictions unchanged.
#' @param known_pi Optional scalar or vector with one known propensity per row.
#'   This skips propensity fitting while retaining cross-fitting for the
#'   marginal outcome regression.
#' @param known_e Optional known `E(A|X)` on the `{-1,+1}` scale. Supply at
#'   most one of `known_pi` and `known_e`.
#' @param verbose Passed to [SuperLearner::SuperLearner()].
#' @param env Environment in which custom Super Learner wrappers are resolved.
#'
#' @return An `odrl_nuisance` object containing aligned out-of-fold
#'   predictions. Its `fits` component records the libraries, coefficients,
#'   actual numbers of inner folds, and warnings for each outer fold. The
#'   inner fold counts are capped by the training sample size or the smaller
#'   treatment arm. The `warnings` component collects distinct warning
#'   messages across folds.
#' @export
odrl_nuisance_sl <- function(
    x, a, y, folds = 5L,
    sl.library = c("SL.mean", "SL.glm", "SL.glmnet"),
    sl.library.pi = sl.library, sl.library.m = sl.library,
    inner_folds = 5L, seed = 1L, positive = NULL,
    propensity_bounds = NULL,
    known_pi = NULL, known_e = NULL,
    verbose = FALSE, env = parent.frame(), fold_id = NULL) {
  .odrl_require("SuperLearner", "for cross-fitted nuisance estimation")
  env <- .odrl_superlearner_env(env)
  encoded <- .odrl_encode_x_fit(x)
  x_matrix <- encoded$x
  treatment_info <- .odrl_encode_treatment(a, positive = positive)
  a_pm1 <- treatment_info$pm1
  treatment_01 <- as.integer(a_pm1 == 1)
  y <- .odrl_validate_outcome(y, nrow(x_matrix))
  if (length(a_pm1) != nrow(x_matrix)) {
    .odrl_abort("`x`, `a`, and `y` must have matching finite rows.")
  }
  if (is.null(fold_id)) {
    .odrl_check_scalar(folds, "folds", 2, nrow(x_matrix), integer = TRUE)
  }
  .odrl_check_scalar(inner_folds, "inner_folds", 2, nrow(x_matrix),
                     integer = TRUE)
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  known_probability <- .odrl_known_propensity(
    known_pi, known_e, nrow(x_matrix)
  )
  estimate_propensity <- is.null(known_probability)

  fold_spec <- .odrl_crossfit_fold_spec(
    a_pm1, folds, fold_id, seed,
    require_training_arms = estimate_propensity
  )
  fold_id <- fold_spec$id
  folds <- fold_spec$folds
  pi_hat <- m_hat <- rep(NA_real_, nrow(x_matrix))
  fit_audit <- vector("list", folds)
  x_frame <- as.data.frame(x_matrix, check.names = FALSE)

  # Fill predictions in their original row positions. Each row is predicted
  # exactly once, from models trained on the other folds.
  for (fold in seq_len(folds)) {
    training_rows <- fold_id != fold
    held_out_rows <- !training_rows
    training_x <- x_frame[training_rows, , drop = FALSE]
    held_out_x <- x_frame[held_out_rows, , drop = FALSE]
    outcome_inner_folds <- min(as.integer(inner_folds), sum(training_rows))

    # Fit P(A = +1 | X) with a binomial model unless the probability is known.
    # Its inner folds are limited by the smaller treatment arm.
    if (estimate_propensity) {
      arm_counts <- table(factor(treatment_01[training_rows], levels = 0:1))
      if (min(arm_counts) < 2L) {
        .odrl_abort(
          "An outer training fold has fewer than two observations in one ",
          "treatment arm; reduce `folds` or use a design with more arm support."
        )
      }
      propensity_inner_folds <- min(
        as.integer(inner_folds), as.integer(min(arm_counts))
      )
      set.seed(seed + 1000L + fold)
      propensity_result <- .odrl_capture_warnings(SuperLearner::SuperLearner(
        Y = treatment_01[training_rows],
        X = training_x,
        newX = held_out_x,
        family = stats::binomial(),
        SL.library = sl.library.pi,
        method = "method.NNloglik",
        cvControl = list(
          V = propensity_inner_folds, stratifyCV = TRUE, shuffle = TRUE
        ),
        control = list(saveFitLibrary = FALSE),
        verbose = verbose,
        env = env
      ))
      propensity_fit <- propensity_result$value
      pi_hat[held_out_rows] <- as.numeric(propensity_fit$SL.predict)
    } else {
      propensity_fit <- NULL
      propensity_result <- list(value = NULL, warnings = character())
      propensity_inner_folds <- NA_integer_
      pi_hat[held_out_rows] <- known_probability[held_out_rows]
    }

    # Fit the pooled mean E(Y | X), using covariates from both treatment arms.
    set.seed(seed + 2000L + fold)
    outcome_result <- .odrl_capture_warnings(SuperLearner::SuperLearner(
      Y = y[training_rows],
      X = training_x,
      newX = held_out_x,
      family = stats::gaussian(),
      SL.library = sl.library.m,
      method = "method.NNLS",
      cvControl = list(V = outcome_inner_folds, shuffle = TRUE),
      control = list(saveFitLibrary = FALSE),
      verbose = verbose,
      env = env
    ))
    outcome_fit <- outcome_result$value
    m_hat[held_out_rows] <- as.numeric(outcome_fit$SL.predict)

    # Retain the model choices and warnings needed to inspect each fold.
    fit_audit[[fold]] <- list(
      fold = fold,
      fold_label = fold_spec$labels[[fold]],
      holdout_n = sum(held_out_rows),
      propensity_library = if (estimate_propensity) {
        propensity_fit$libraryNames
      } else {
        "known by design"
      },
      propensity_coef = if (estimate_propensity) {
        propensity_fit$coef
      } else {
        NA_real_
      },
      outcome_library = outcome_fit$libraryNames,
      outcome_coef = outcome_fit$coef,
      propensity_inner_folds = propensity_inner_folds,
      outcome_inner_folds = outcome_inner_folds,
      propensity_warnings = propensity_result$warnings,
      outcome_warnings = outcome_result$warnings
    )
  }
  if (any(!is.finite(c(pi_hat, m_hat)))) {
    .odrl_abort("Super Learner returned non-finite nuisance predictions.")
  }

  # Apply any requested probability bounds after collecting all predictions.
  nuisance <- odrl_nuisance_user(
    m = m_hat,
    pi = pi_hat,
    fold_id = fold_id,
    source = paste0(
      folds, "-fold cross-fitted marginal outcome Super Learner; propensity ",
      if (estimate_propensity) {
        paste0("Super Learner: ", .odrl_sl_library_label(sl.library.pi))
      } else {
        "known by design"
      },
      "; outcome library: ", .odrl_sl_library_label(sl.library.m)
    ),
    out_of_fold = TRUE,
    propensity_bounds = propensity_bounds
  )
  nuisance$fits <- fit_audit
  nuisance$blueprint <- encoded$blueprint
  nuisance$treatment_map <- treatment_info
  nuisance$sl.library <- sl.library
  nuisance$sl.library.pi <- sl.library.pi
  nuisance$sl.library.m <- sl.library.m
  nuisance$known_propensity <- !estimate_propensity
  nuisance$folds <- folds
  nuisance$fold_source <- fold_spec$source
  nuisance$fold_labels <- fold_spec$labels
  nuisance$inner_folds <- inner_folds
  nuisance$warnings <- unique(unlist(lapply(fit_audit, function(fit) {
    c(fit$propensity_warnings, fit$outcome_warnings)
  }), use.names = FALSE))
  nuisance
}

#' Estimate treatment and outcome predictions with regression
#'
#' Fits logistic regression for the treatment probability and linear
#' regression for the mean outcome. Both models include an intercept and the
#' columns of `x`. The outcome model pools both treatment arms and uses
#' covariates alone. Supply additional columns in `x` to include terms such
#' as squares or interactions.
#'
#' Each fold is held out in turn while models are fitted on the remaining
#' rows. Predictions are returned in the original row order and can be passed
#' directly to [odrl()].
#'
#' @inheritParams odrl_nuisance_sl
#'
#' @return An `odrl_nuisance` object containing aligned out-of-fold
#'   predictions and coefficient, rank, convergence, and warning diagnostics
#'   for each fold.
#' @export
odrl_nuisance_parametric <- function(
    x, a, y, folds = 5L, fold_id = NULL, seed = 1L, positive = NULL,
    propensity_bounds = NULL, known_pi = NULL, known_e = NULL) {
  encoded <- .odrl_encode_x_fit(x)
  x_matrix <- encoded$x
  treatment_info <- .odrl_encode_treatment(a, positive = positive)
  a_pm1 <- treatment_info$pm1
  treatment_01 <- as.integer(a_pm1 == 1)
  y <- .odrl_validate_outcome(y, nrow(x_matrix))
  if (length(a_pm1) != nrow(x_matrix)) {
    .odrl_abort("`x`, `a`, and `y` must have matching finite rows.")
  }
  if (is.null(fold_id)) {
    .odrl_check_scalar(folds, "folds", 2, nrow(x_matrix), integer = TRUE)
  }
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  known_probability <- .odrl_known_propensity(
    known_pi, known_e, nrow(x_matrix)
  )
  estimate_propensity <- is.null(known_probability)

  fold_spec <- .odrl_crossfit_fold_spec(
    a_pm1, folds, fold_id, seed,
    require_training_arms = estimate_propensity
  )
  fold_id <- fold_spec$id
  folds <- fold_spec$folds
  design <- cbind(`(Intercept)` = 1, x_matrix)
  pi_hat <- m_hat <- rep(NA_real_, nrow(x_matrix))
  fit_audit <- vector("list", folds)

  for (fold in seq_len(folds)) {
    training_rows <- fold_id != fold
    held_out_rows <- !training_rows
    training_design <- design[training_rows, , drop = FALSE]
    held_out_design <- design[held_out_rows, , drop = FALSE]

    # Logistic regression predicts the probability of the positive treatment.
    if (estimate_propensity) {
      propensity_result <- .odrl_capture_warnings(stats::glm.fit(
        x = training_design, y = treatment_01[training_rows],
        family = stats::binomial()
      ))
      propensity_fit <- propensity_result$value
      propensity_coef <- propensity_fit$coefficients

      # glm.fit marks redundant columns with NA coefficients. Giving those
      # columns zero contribution preserves the fitted linear predictor.
      propensity_aliased <- is.na(propensity_coef)
      propensity_coef[propensity_aliased] <- 0
      pi_hat[held_out_rows] <- stats::plogis(drop(
        held_out_design %*% propensity_coef
      ))
    } else {
      propensity_fit <- NULL
      propensity_result <- list(value = NULL, warnings = character())
      propensity_coef <- rep(NA_real_, ncol(design))
      names(propensity_coef) <- colnames(design)
      propensity_aliased <- rep(FALSE, ncol(design))
      pi_hat[held_out_rows] <- known_probability[held_out_rows]
    }

    # Linear regression predicts E(Y | X), pooling both treatment arms.
    outcome_result <- .odrl_capture_warnings(stats::lm.fit(
      x = training_design, y = y[training_rows]
    ))
    outcome_fit <- outcome_result$value
    outcome_coef <- outcome_fit$coefficients
    outcome_aliased <- is.na(outcome_coef)
    outcome_coef[outcome_aliased] <- 0
    m_hat[held_out_rows] <- drop(held_out_design %*% outcome_coef)

    fit_audit[[fold]] <- list(
      fold = fold,
      fold_label = fold_spec$labels[[fold]],
      holdout_n = sum(held_out_rows),
      propensity_engine = if (estimate_propensity) {
        "binomial logistic regression"
      } else {
        "known by design"
      },
      propensity_coef = propensity_coef,
      propensity_rank = if (estimate_propensity) {
        propensity_fit$rank
      } else {
        NA_integer_
      },
      propensity_converged = if (estimate_propensity) {
        propensity_fit$converged
      } else {
        NA
      },
      propensity_aliased = names(propensity_coef)[propensity_aliased],
      outcome_engine = "Gaussian linear regression",
      outcome_coef = outcome_coef,
      outcome_rank = outcome_fit$rank,
      outcome_aliased = names(outcome_coef)[outcome_aliased],
      propensity_warnings = propensity_result$warnings,
      outcome_warnings = outcome_result$warnings
    )
  }
  if (any(!is.finite(c(pi_hat, m_hat)))) {
    .odrl_abort("Parametric models returned non-finite nuisance predictions.")
  }

  nuisance <- odrl_nuisance_user(
    m = m_hat,
    pi = pi_hat,
    fold_id = fold_id,
    source = paste0(
      folds, "-fold cross-fitted main-effects linear outcome regression; ",
      "propensity ", if (estimate_propensity) {
        "main-effects logistic regression"
      } else {
        "known by design"
      }
    ),
    out_of_fold = TRUE,
    propensity_bounds = propensity_bounds
  )
  nuisance$fits <- fit_audit
  nuisance$blueprint <- encoded$blueprint
  nuisance$treatment_map <- treatment_info
  nuisance$engine <- "parametric"
  nuisance$known_propensity <- !estimate_propensity
  nuisance$folds <- folds
  nuisance$fold_source <- fold_spec$source
  nuisance$fold_labels <- fold_spec$labels
  nuisance$warnings <- unique(unlist(lapply(fit_audit, function(fit) {
    c(fit$propensity_warnings, fit$outcome_warnings)
  }), use.names = FALSE))
  nuisance
}

.odrl_capture_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(messages))
}

.odrl_resolve_nuisance <- function(
    nuisance, x, a, y, nuisance_folds, sl.library, sl.library.pi,
    sl.library.m, sl_inner_folds, propensity_bounds, known_pi, known_e,
    seed, positive, sl_verbose, sl_env, nuisance_fold_id = NULL) {
  built_in_sl <- is.null(nuisance) || identical(nuisance, "superlearner")
  built_in_parametric <- identical(nuisance, "parametric") ||
    identical(nuisance, "glm")
  built_in <- built_in_sl || built_in_parametric

  # Known propensities and requested folds configure model fitting. Supplied
  # predictions already carry their own values and fold information.
  if (!built_in && (!is.null(known_pi) || !is.null(known_e))) {
    .odrl_abort(
      "`known_pi` and `known_e` are only used with built-in nuisance ",
      "estimation. Do not combine them with a supplied nuisance object."
    )
  }
  if (!built_in && !is.null(nuisance_fold_id)) {
    .odrl_abort(
      "`nuisance_fold_id` is only used with built-in nuisance estimation. ",
      "Construct supplied nuisance predictions with their own fold metadata."
    )
  }
  if (built_in_sl) {
    return(odrl_nuisance_sl(
      x = x, a = a, y = y, folds = nuisance_folds,
      fold_id = nuisance_fold_id,
      sl.library = sl.library, sl.library.pi = sl.library.pi,
      sl.library.m = sl.library.m, inner_folds = sl_inner_folds,
      seed = seed, positive = positive,
      propensity_bounds = propensity_bounds,
      known_pi = known_pi, known_e = known_e,
      verbose = sl_verbose, env = sl_env
    ))
  }
  if (built_in_parametric) {
    return(odrl_nuisance_parametric(
      x = x, a = a, y = y, folds = nuisance_folds,
      fold_id = nuisance_fold_id,
      seed = seed, positive = positive,
      propensity_bounds = propensity_bounds,
      known_pi = known_pi, known_e = known_e
    ))
  }
  if (inherits(nuisance, "odrl_nuisance")) {
    if (!is.null(propensity_bounds) &&
        !identical(propensity_bounds, nuisance$propensity_bounds)) {
      .odrl_abort(
        "`propensity_bounds` cannot silently override a constructed nuisance ",
        "object. Apply the bounds in `odrl_nuisance_user()` or ",
        "`odrl_nuisance_sl()` before fitting."
      )
    }
    return(nuisance)
  }
  if (is.function(nuisance)) {
    supplied <- nuisance(x = x, a = a, y = y)
    return(.odrl_resolve_nuisance(
      nuisance = supplied, x = x, a = a, y = y,
      nuisance_folds = nuisance_folds,
      sl.library = sl.library, sl.library.pi = sl.library.pi,
      sl.library.m = sl.library.m, sl_inner_folds = sl_inner_folds,
      propensity_bounds = propensity_bounds,
      known_pi = known_pi, known_e = known_e,
      seed = seed, positive = positive,
      sl_verbose = sl_verbose, sl_env = sl_env,
      nuisance_fold_id = nuisance_fold_id
    ))
  }
  if (is.list(nuisance)) {
    return(odrl_nuisance_user(
      m = nuisance$m,
      pi = nuisance$pi,
      e = nuisance$e,
      fold_id = nuisance$fold_id,
      source = nuisance$source %||% "user supplied list",
      out_of_fold = nuisance$out_of_fold %||% FALSE,
      propensity_bounds = propensity_bounds
    ))
  }
  .odrl_abort(
    "`nuisance` must be NULL, \"superlearner\", \"parametric\", an ",
    "`odrl_nuisance` ",
    "object, or a list containing `m` and exactly one of `pi` or `e`."
  )
}

.odrl_sl_library_label <- function(library) {
  entry <- if (is.character(library)) as.list(library) else library
  paste(vapply(entry, function(x) paste(as.character(x), collapse = "/"),
               character(1)), collapse = ", ")
}

.odrl_superlearner_env <- function(user_env) {
  if (!is.environment(user_env)) {
    .odrl_abort("`env` must be an environment.")
  }
  resolver <- new.env(parent = user_env)
  exports <- getNamespaceExports("SuperLearner")
  for (name in exports) {
    if (!exists(name, envir = resolver, inherits = TRUE)) {
      assign(name, getExportedValue("SuperLearner", name), envir = resolver)
    }
  }
  resolver
}

`%||%` <- function(x, y) if (is.null(x)) y else x
