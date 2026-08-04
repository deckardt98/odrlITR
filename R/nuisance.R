#' Construct user-supplied ODRL nuisance predictions
#'
#' User-supplied predictions should be aligned with the rows of `x`. Unless a
#' nuisance is known by design, they should be out-of-fold predictions. Supply
#' either `pi` or `e`, where `e = 2 * pi - 1` under the package's internal
#' \eqn{\{-1,+1\}} treatment coding.
#'
#' @param m Numeric predictions of \eqn{E(Y\mid X)}.
#' @param pi Numeric predictions of \eqn{P(A=+1\mid X)}.
#' @param e Numeric predictions of \eqn{E(A\mid X)}.
#' @param fold_id Optional outer-fold identifier for auditing.
#' @param source Short description of how the predictions were obtained.
#' @param out_of_fold Whether predictions are out of fold or known by design.
#'   The conservative default is `FALSE`; [odrl()] warns because ordinary
#'   in-sample predictions do not implement the paper's cross-fitting step.
#' @param propensity_bounds Optional two-element numerical safeguard. Use
#'   `NULL` to retain supplied probabilities exactly.
#'
#' @return An object of class `odrl_nuisance`.
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
  label <- as.character(fold_id)
  if (any(!nzchar(trimws(label)))) {
    .odrl_abort("`fold_id` cannot contain empty labels.")
  }
  levels <- unique(label)
  if (length(levels) < 2L) {
    .odrl_abort("`fold_id` must identify at least two nonempty folds.")
  }
  id <- match(label, levels)
  if (min_training_n > 0L) {
    training_sizes <- vapply(
      seq_along(levels), function(fold) sum(id != fold), integer(1)
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
    for (fold in seq_along(levels)) {
      counts <- table(factor(a_pm1[id != fold], levels = c(-1, 1)))
      if (min(counts) < 2L) {
        .odrl_abort(
          "Each user-specified outer training set must contain at least two ",
          "observations from each treatment arm."
        )
      }
    }
  }
  list(
    id = if (normalize) as.integer(id) else fold_id,
    labels = levels,
    folds = length(levels),
    source = "user supplied"
  )
}

.odrl_treatment_stratified_folds <- function(
    a_pm1, folds, seed, require_training_arms = TRUE) {
  folds <- as.integer(folds)
  .odrl_check_scalar(folds, "folds", 2, length(a_pm1), integer = TRUE)
  counts <- table(factor(a_pm1, levels = c(-1, 1)))
  if (require_training_arms && min(counts) < folds) {
    .odrl_abort("Each treatment arm must contain at least `folds` rows.")
  }

  set.seed(seed)
  assignment <- integer(length(a_pm1))
  total_load <- integer(folds)
  for (value in c(-1, 1)) {
    index <- sample(which(a_pm1 == value))
    base <- length(index) %/% folds
    remainder <- length(index) %% folds
    arm_load <- rep.int(base, folds)
    if (remainder > 0L) {
      # Put this arm's remainder into the currently smallest folds. Random
      # tie-breaking preserves exchangeability while balancing total size.
      priority <- order(total_load, stats::runif(folds))
      arm_load[priority[seq_len(remainder)]] <-
        arm_load[priority[seq_len(remainder)]] + 1L
    }
    labels <- rep.int(seq_len(folds), arm_load)
    assignment[index] <- sample(labels, length(labels), replace = FALSE)
    total_load <- total_load + arm_load
  }
  assignment
}

.odrl_crossfit_fold_spec <- function(a_pm1, folds, fold_id, seed,
                                     require_training_arms = TRUE) {
  if (is.null(fold_id)) {
    id <- .odrl_treatment_stratified_folds(
      a_pm1, folds, seed,
      require_training_arms = require_training_arms
    )
    if (any(vapply(
      seq_len(as.integer(folds)), function(fold) sum(id != fold), integer(1)
    ) < 2L)) {
      .odrl_abort(
        "Each outer training set must contain at least two observations; ",
        "reduce `folds`."
      )
    }
    if (require_training_arms) {
      for (fold in seq_len(as.integer(folds))) {
        counts <- table(factor(a_pm1[id != fold], levels = c(-1, 1)))
        if (min(counts) < 2L) {
          .odrl_abort(
            "Each outer training set must contain at least two observations ",
            "from each treatment arm; reduce `folds`."
          )
        }
      }
    }
    return(list(
      id = id,
      labels = as.character(seq_len(as.integer(folds))),
      folds = as.integer(folds),
      source = "treatment-stratified"
    ))
  }
  .odrl_validate_fold_id(
    fold_id, length(a_pm1), a_pm1 = a_pm1,
    require_training_arms = require_training_arms,
    min_training_n = 2L, normalize = TRUE
  )
}

#' Cross-fitted Super Learner nuisances for ODRL
#'
#' Fits the propensity \eqn{P(A=+1\mid X)} and pooled marginal outcome
#' regression \eqn{E(Y\mid X)} on outer training folds and predicts held-out
#' rows. Super Learner's internal cross-validation is separate from this outer
#' cross-fitting layer.
#'
#' @param x Covariate matrix or data frame.
#' @param a Binary treatment, coded as `{0,1}`, `{-1,+1}`, logical, or a
#'   two-level factor.
#' @param y Numeric outcome.
#' @param folds Number of outer cross-fitting folds. Ignored when `fold_id` is
#'   supplied. Automatically generated folds are treatment-stratified, and
#'   their total sizes differ by at most one row.
#' @param fold_id Optional row-aligned user-specified outer-fold identifiers.
#'   Labels may be numeric, character, or factor and are normalized internally.
#'   Every resulting training set must contain at least two observations and
#'   must retain at least two observations from each treatment arm when the
#'   propensity is estimated.
#' @param sl.library Character vector or list understood by
#'   [SuperLearner::SuperLearner()].
#' @param sl.library.pi Optional propensity-specific library. Defaults to
#'   `sl.library`.
#' @param sl.library.m Optional marginal-outcome-specific library. Defaults to
#'   `sl.library`.
#' @param inner_folds Super Learner cross-validation folds within each outer
#'   training sample.
#' @param seed Reproducibility seed.
#' @param positive For a factor or character treatment, the level representing
#'   treatment `+1`. By default, the second factor level is positive.
#' @param propensity_bounds Optional numerical safeguard. `NULL` retains the
#'   Super Learner predictions exactly.
#' @param known_pi Optional scalar or row-aligned vector containing a known
#'   propensity. This skips propensity fitting while retaining cross-fitting
#'   for the marginal outcome regression.
#' @param known_e Optional known `E(A|X)` on the `{-1,+1}` scale. Supply at
#'   most one of `known_pi` and `known_e`.
#' @param verbose Passed to [SuperLearner::SuperLearner()].
#' @param env Environment in which custom Super Learner wrappers are resolved.
#'
#' @return An `odrl_nuisance` object containing aligned out-of-fold
#'   predictions. Its `fits` component records the libraries, coefficients,
#'   actual (possibly capped) inner-fold counts, and captured warnings for
#'   every outer fold; `warnings` is the unique aggregate warning ledger.
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
  a_info <- .odrl_encode_treatment(a, positive = positive)
  a_pm1 <- a_info$pm1
  a01 <- as.integer(a_pm1 == 1)
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
  if (!is.null(known_pi) && !is.null(known_e)) {
    .odrl_abort("Supply at most one of `known_pi` and `known_e`.")
  }
  known_probability <- NULL
  if (!is.null(known_pi) || !is.null(known_e)) {
    known_input <- if (!is.null(known_pi)) known_pi else known_e
    if (!is.numeric(known_input) || is.factor(known_input)) {
      .odrl_abort("Known propensity values must be genuinely numeric.")
    }
    known_probability <- if (!is.null(known_pi)) {
      as.numeric(known_pi)
    } else {
      (as.numeric(known_e) + 1) / 2
    }
    if (length(known_probability) == 1L) {
      known_probability <- rep(known_probability, nrow(x_matrix))
    }
    if (length(known_probability) != nrow(x_matrix) ||
        any(!is.finite(known_probability)) ||
        any(known_probability < 0 | known_probability > 1)) {
      .odrl_abort("Known propensity values must align with rows and lie in [0,1].")
    }
  }
  fold_spec <- .odrl_crossfit_fold_spec(
    a_pm1, folds, fold_id, seed,
    require_training_arms = is.null(known_probability)
  )
  fold_id <- fold_spec$id
  folds <- fold_spec$folds
  pi_hat <- m_hat <- rep(NA_real_, nrow(x_matrix))
  fit_audit <- vector("list", folds)
  x_frame <- as.data.frame(x_matrix, check.names = FALSE)
  for (fold in seq_len(folds)) {
    train <- fold_id != fold
    holdout <- !train
    inner_v_m <- min(as.integer(inner_folds), sum(train))
    if (is.null(known_probability)) {
      arm_counts <- table(factor(a01[train], levels = 0:1))
      if (min(arm_counts) < 2L) {
        .odrl_abort(
          "An outer training fold has fewer than two observations in one ",
          "treatment arm; reduce `folds` or use a design with more arm support."
        )
      }
      inner_v_pi <- min(as.integer(inner_folds),
                        as.integer(min(arm_counts)))
      set.seed(seed + 1000L + fold)
      pi_result <- .odrl_capture_warnings(SuperLearner::SuperLearner(
        Y = a01[train],
        X = x_frame[train, , drop = FALSE],
        newX = x_frame[holdout, , drop = FALSE],
        family = stats::binomial(),
        SL.library = sl.library.pi,
        method = "method.NNloglik",
        cvControl = list(V = inner_v_pi, stratifyCV = TRUE, shuffle = TRUE),
        control = list(saveFitLibrary = FALSE),
        verbose = verbose,
        env = env
      ))
      pi_fit <- pi_result$value
      pi_hat[holdout] <- as.numeric(pi_fit$SL.predict)
    } else {
      pi_fit <- NULL
      pi_result <- list(value = NULL, warnings = character())
      pi_hat[holdout] <- known_probability[holdout]
    }
    set.seed(seed + 2000L + fold)
    m_result <- .odrl_capture_warnings(SuperLearner::SuperLearner(
      Y = y[train],
      X = x_frame[train, , drop = FALSE],
      newX = x_frame[holdout, , drop = FALSE],
      family = stats::gaussian(),
      SL.library = sl.library.m,
      method = "method.NNLS",
      cvControl = list(V = inner_v_m, shuffle = TRUE),
      control = list(saveFitLibrary = FALSE),
      verbose = verbose,
      env = env
    ))
    m_fit <- m_result$value
    m_hat[holdout] <- as.numeric(m_fit$SL.predict)
    fit_audit[[fold]] <- list(
      fold = fold,
      fold_label = fold_spec$labels[[fold]],
      holdout_n = sum(holdout),
      propensity_library = if (is.null(pi_fit)) "known by design" else {
        pi_fit$libraryNames
      },
      propensity_coef = if (is.null(pi_fit)) NA_real_ else pi_fit$coef,
      outcome_library = m_fit$libraryNames,
      outcome_coef = m_fit$coef,
      propensity_inner_folds = if (is.null(pi_fit)) NA_integer_ else inner_v_pi,
      outcome_inner_folds = inner_v_m,
      propensity_warnings = pi_result$warnings,
      outcome_warnings = m_result$warnings
    )
  }
  if (any(!is.finite(c(pi_hat, m_hat)))) {
    .odrl_abort("Super Learner returned non-finite nuisance predictions.")
  }
  nuisance <- odrl_nuisance_user(
    m = m_hat,
    pi = pi_hat,
    fold_id = fold_id,
    source = paste0(
      folds, "-fold cross-fitted marginal outcome Super Learner; propensity ",
      if (is.null(known_probability)) {
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
  nuisance$treatment_map <- a_info
  nuisance$sl.library <- sl.library
  nuisance$sl.library.pi <- sl.library.pi
  nuisance$sl.library.m <- sl.library.m
  nuisance$known_propensity <- !is.null(known_probability)
  nuisance$folds <- folds
  nuisance$fold_source <- fold_spec$source
  nuisance$fold_labels <- fold_spec$labels
  nuisance$inner_folds <- inner_folds
  nuisance$warnings <- unique(unlist(lapply(fit_audit, function(fit) {
    c(fit$propensity_warnings, fit$outcome_warnings)
  }), use.names = FALSE))
  nuisance
}

#' Cross-fitted parametric nuisances for ODRL
#'
#' Fits a main-effects logistic regression for the propensity and a
#' main-effects linear regression for the pooled marginal outcome. Models are
#' trained separately within each outer training set and predict only held-out
#' observations. To use a richer prespecified basis, pass that basis as `x`.
#'
#' @inheritParams odrl_nuisance_sl
#'
#' @return An `odrl_nuisance` object containing aligned out-of-fold
#'   predictions and fold-specific coefficient, rank, convergence, and warning
#'   diagnostics.
#' @export
odrl_nuisance_parametric <- function(
    x, a, y, folds = 5L, fold_id = NULL, seed = 1L, positive = NULL,
    propensity_bounds = NULL, known_pi = NULL, known_e = NULL) {
  encoded <- .odrl_encode_x_fit(x)
  x_matrix <- encoded$x
  a_info <- .odrl_encode_treatment(a, positive = positive)
  a_pm1 <- a_info$pm1
  a01 <- as.integer(a_pm1 == 1)
  y <- .odrl_validate_outcome(y, nrow(x_matrix))
  if (length(a_pm1) != nrow(x_matrix)) {
    .odrl_abort("`x`, `a`, and `y` must have matching finite rows.")
  }
  if (is.null(fold_id)) {
    .odrl_check_scalar(folds, "folds", 2, nrow(x_matrix), integer = TRUE)
  }
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  if (!is.null(known_pi) && !is.null(known_e)) {
    .odrl_abort("Supply at most one of `known_pi` and `known_e`.")
  }
  known_probability <- NULL
  if (!is.null(known_pi) || !is.null(known_e)) {
    known_input <- if (!is.null(known_pi)) known_pi else known_e
    if (!is.numeric(known_input) || is.factor(known_input)) {
      .odrl_abort("Known propensity values must be genuinely numeric.")
    }
    known_probability <- if (!is.null(known_pi)) {
      as.numeric(known_pi)
    } else {
      (as.numeric(known_e) + 1) / 2
    }
    if (length(known_probability) == 1L) {
      known_probability <- rep(known_probability, nrow(x_matrix))
    }
    if (length(known_probability) != nrow(x_matrix) ||
        any(!is.finite(known_probability)) ||
        any(known_probability < 0 | known_probability > 1)) {
      .odrl_abort(
        "Known propensity values must align with rows and lie in [0,1]."
      )
    }
  }

  fold_spec <- .odrl_crossfit_fold_spec(
    a_pm1, folds, fold_id, seed,
    require_training_arms = is.null(known_probability)
  )
  fold_id <- fold_spec$id
  folds <- fold_spec$folds
  design <- cbind(`(Intercept)` = 1, x_matrix)
  pi_hat <- m_hat <- rep(NA_real_, nrow(x_matrix))
  fit_audit <- vector("list", folds)

  for (fold in seq_len(folds)) {
    train <- fold_id != fold
    holdout <- !train
    if (is.null(known_probability)) {
      pi_result <- .odrl_capture_warnings(stats::glm.fit(
        x = design[train, , drop = FALSE], y = a01[train],
        family = stats::binomial()
      ))
      pi_fit <- pi_result$value
      pi_coef <- pi_fit$coefficients
      pi_aliased <- is.na(pi_coef)
      pi_coef[pi_aliased] <- 0
      pi_hat[holdout] <- stats::plogis(drop(
        design[holdout, , drop = FALSE] %*% pi_coef
      ))
    } else {
      pi_fit <- NULL
      pi_result <- list(value = NULL, warnings = character())
      pi_coef <- rep(NA_real_, ncol(design))
      names(pi_coef) <- colnames(design)
      pi_aliased <- rep(FALSE, ncol(design))
      pi_hat[holdout] <- known_probability[holdout]
    }

    m_result <- .odrl_capture_warnings(stats::lm.fit(
      x = design[train, , drop = FALSE], y = y[train]
    ))
    m_fit <- m_result$value
    m_coef <- m_fit$coefficients
    m_aliased <- is.na(m_coef)
    m_coef[m_aliased] <- 0
    m_hat[holdout] <- drop(
      design[holdout, , drop = FALSE] %*% m_coef
    )

    fit_audit[[fold]] <- list(
      fold = fold,
      fold_label = fold_spec$labels[[fold]],
      holdout_n = sum(holdout),
      propensity_engine = if (is.null(pi_fit)) {
        "known by design"
      } else {
        "binomial logistic regression"
      },
      propensity_coef = pi_coef,
      propensity_rank = if (is.null(pi_fit)) NA_integer_ else pi_fit$rank,
      propensity_converged = if (is.null(pi_fit)) NA else pi_fit$converged,
      propensity_aliased = names(pi_coef)[pi_aliased],
      outcome_engine = "Gaussian linear regression",
      outcome_coef = m_coef,
      outcome_rank = m_fit$rank,
      outcome_aliased = names(m_coef)[m_aliased],
      propensity_warnings = pi_result$warnings,
      outcome_warnings = m_result$warnings
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
      "propensity ", if (is.null(known_probability)) {
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
  nuisance$treatment_map <- a_info
  nuisance$engine <- "parametric"
  nuisance$known_propensity <- !is.null(known_probability)
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
      supplied, x, a, y, nuisance_folds, sl.library, sl.library.pi,
      sl.library.m, sl_inner_folds, propensity_bounds, known_pi, known_e,
      seed, positive, sl_verbose, sl_env, nuisance_fold_id
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
  namespace <- asNamespace("SuperLearner")
  exports <- getNamespaceExports("SuperLearner")
  for (name in exports) {
    if (!exists(name, envir = resolver, inherits = TRUE)) {
      assign(name, getExportedValue("SuperLearner", name), envir = resolver)
    }
  }
  resolver
}

`%||%` <- function(x, y) if (is.null(x)) y else x
