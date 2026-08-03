.odrl_abort <- function(...) {
  stop(..., call. = FALSE)
}

.odrl_require <- function(package, reason) {
  if (!requireNamespace(package, quietly = TRUE)) {
    .odrl_abort(
      "Package '", package, "' is required ", reason,
      ". Install it with install.packages(\"", package, "\")."
    )
  }
}

.odrl_validate_outcome <- function(y, n = NULL) {
  if (!is.numeric(y) || is.factor(y)) {
    .odrl_abort("`y` must be a numeric outcome, not factor or character data.")
  }
  y <- as.numeric(y)
  if ((!is.null(n) && length(y) != n) || !length(y) || any(!is.finite(y))) {
    .odrl_abort("`y` must contain one finite numeric value per row of `x`.")
  }
  y
}

.odrl_check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                               integer = FALSE, open_lower = FALSE,
                               open_upper = FALSE) {
  if (length(x) != 1L || !is.finite(x)) {
    .odrl_abort("`", name, "` must be one finite number.")
  }
  if (integer && x != as.integer(x)) {
    .odrl_abort("`", name, "` must be an integer.")
  }
  lower_bad <- if (open_lower) x <= lower else x < lower
  upper_bad <- if (open_upper) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    interval <- paste0(
      if (open_lower) "(" else "[", lower, ", ", upper,
      if (open_upper) ")" else "]"
    )
    .odrl_abort("`", name, "` must lie in ", interval, ".")
  }
  invisible(TRUE)
}

.odrl_encode_treatment <- function(a, positive = NULL) {
  if (anyNA(a)) .odrl_abort("`a` cannot contain missing values.")
  if (is.logical(a)) {
    return(list(
      pm1 = ifelse(a, 1, -1),
      positive = TRUE,
      negative = FALSE,
      original = "logical"
    ))
  }
  if (is.factor(a) || is.character(a)) {
    f <- if (is.factor(a)) droplevels(a) else factor(a)
    if (nlevels(f) != 2L) .odrl_abort("`a` must have exactly two levels.")
    lev <- levels(f)
    if (!is.null(positive)) {
      if (length(positive) != 1L || is.na(positive) ||
          !nzchar(as.character(positive))) {
        .odrl_abort("`positive` must be one nonempty treatment level.")
      }
      positive <- as.character(positive)
      if (!positive %in% lev) {
        .odrl_abort("`positive` must be one of the two treatment levels.")
      }
      lev <- c(setdiff(lev, positive), positive)
      f <- factor(as.character(f), levels = lev)
    }
    return(list(
      pm1 = ifelse(f == lev[[2L]], 1, -1),
      positive = lev[[2L]],
      negative = lev[[1L]],
      original = if (is.factor(a)) "factor" else "character"
    ))
  }
  if (!is.numeric(a) && !is.integer(a)) {
    .odrl_abort("`a` must be binary numeric, logical, character, or factor.")
  }
  values <- sort(unique(as.numeric(a)))
  if (identical(values, c(-1, 1))) {
    return(list(pm1 = as.numeric(a), positive = 1, negative = -1,
                original = "pm1"))
  }
  if (identical(values, c(0, 1))) {
    return(list(pm1 = 2 * as.numeric(a) - 1, positive = 1, negative = 0,
                original = "zero_one"))
  }
  .odrl_abort("Numeric `a` must contain exactly {0,1} or {-1,+1}.")
}

.odrl_decode_action <- function(action, treatment_map) {
  positive <- treatment_map$positive
  negative <- treatment_map$negative
  if (identical(treatment_map$original, "factor")) {
    return(factor(
      ifelse(action == 1, positive, negative),
      levels = c(negative, positive)
    ))
  }
  if (identical(treatment_map$original, "character")) {
    return(ifelse(action == 1, positive, negative))
  }
  if (identical(treatment_map$original, "logical")) {
    return(action == 1)
  }
  if (identical(treatment_map$original, "pm1")) return(as.numeric(action))
  if (identical(treatment_map$original, "zero_one")) {
    return(as.numeric(action == 1))
  }
  ifelse(action == 1, positive, negative)
}

.odrl_encode_x_fit <- function(x) {
  if (is.data.frame(x)) {
    if (!ncol(x)) .odrl_abort("`x` must contain at least one covariate.")
    formula <- stats::as.formula("~ . - 1")
    terms <- stats::terms(formula, data = x)
    matrix <- stats::model.matrix(terms, data = x, na.action = stats::na.pass)
    blueprint <- list(
      kind = "data.frame",
      terms = terms,
      xlevels = stats::.getXlevels(terms, x),
      contrasts = attr(matrix, "contrasts"),
      columns = colnames(matrix)
    )
  } else {
    matrix <- as.matrix(x)
    if (is.null(dim(matrix)) || !ncol(matrix)) {
      .odrl_abort("`x` must be a matrix or data frame with covariates.")
    }
    storage.mode(matrix) <- "double"
    if (is.null(colnames(matrix))) {
      colnames(matrix) <- paste0("X", seq_len(ncol(matrix)))
    }
    blueprint <- list(kind = "matrix", columns = colnames(matrix))
  }
  storage.mode(matrix) <- "double"
  if (!nrow(matrix) || any(!is.finite(matrix))) {
    .odrl_abort("The encoded covariate matrix must be nonempty and finite.")
  }
  if (anyDuplicated(colnames(matrix))) {
    .odrl_abort("Encoded covariate names must be unique.")
  }
  list(x = matrix, blueprint = blueprint)
}

.odrl_encode_x_new <- function(x, blueprint) {
  if (identical(blueprint$kind, "formula")) {
    if (!is.data.frame(x)) {
      .odrl_abort("Formula fits require `newdata` to be a data frame.")
    }
    if (!is.null(blueprint$treatment) && blueprint$treatment %in% names(x)) {
      x[[blueprint$treatment]] <- NULL
    }
    frame <- stats::model.frame(
      blueprint$terms, data = x, xlev = blueprint$xlevels,
      na.action = stats::na.fail
    )
    matrix <- stats::model.matrix(
      blueprint$terms, data = frame,
      contrasts.arg = blueprint$contrasts
    )
    matrix <- matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
    missing <- setdiff(blueprint$columns, colnames(matrix))
    extra <- setdiff(colnames(matrix), blueprint$columns)
    if (length(missing) || length(extra)) {
      .odrl_abort("`newdata` does not match the fitted formula encoding.")
    }
    matrix <- matrix[, blueprint$columns, drop = FALSE]
  } else if (identical(blueprint$kind, "data.frame")) {
    if (!is.data.frame(x)) x <- as.data.frame(x)
    matrix <- stats::model.matrix(
      blueprint$terms,
      data = x,
      contrasts.arg = blueprint$contrasts,
      xlev = blueprint$xlevels,
      na.action = stats::na.pass
    )
    missing <- setdiff(blueprint$columns, colnames(matrix))
    extra <- setdiff(colnames(matrix), blueprint$columns)
    if (length(missing) || length(extra)) {
      .odrl_abort("`newdata` does not match the training covariate encoding.")
    }
    matrix <- matrix[, blueprint$columns, drop = FALSE]
  } else {
    matrix <- as.matrix(x)
    storage.mode(matrix) <- "double"
    if (ncol(matrix) != length(blueprint$columns)) {
      .odrl_abort("`newdata` has a different number of covariates.")
    }
    if (!is.null(colnames(matrix))) {
      if (!setequal(colnames(matrix), blueprint$columns)) {
        .odrl_abort("Named `newdata` columns do not match training columns.")
      }
      matrix <- matrix[, blueprint$columns, drop = FALSE]
    }
  }
  storage.mode(matrix) <- "double"
  if (any(!is.finite(matrix))) .odrl_abort("`newdata` must be finite.")
  matrix
}

.odrl_standardize_fit <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-10] <- 1
  list(center = center, scale = scale)
}

.odrl_standardize_apply <- function(x, transform) {
  sweep(sweep(as.matrix(x), 2L, transform$center, "-"),
        2L, transform$scale, "/")
}

.odrl_minmax_fit <- function(x) {
  lower <- apply(x, 2L, min)
  upper <- apply(x, 2L, max)
  span <- upper - lower
  span[!is.finite(span) | span < 1e-10] <- 1
  list(lower = lower, span = span)
}

.odrl_minmax_apply <- function(x, transform) {
  sweep(sweep(as.matrix(x), 2L, transform$lower, "-"),
        2L, transform$span, "/")
}

.odrl_balanced_folds <- function(a_pm1, folds, seed) {
  folds <- as.integer(folds)
  .odrl_check_scalar(folds, "folds", 2, length(a_pm1), integer = TRUE)
  assignment <- integer(length(a_pm1))
  set.seed(seed)
  for (value in c(-1, 1)) {
    index <- sample(which(a_pm1 == value))
    if (length(index) < folds) {
      .odrl_abort("Each treatment arm must contain at least `folds` rows.")
    }
    assignment[index] <- rep(seq_len(folds), length.out = length(index))
  }
  assignment
}

.odrl_score_folds <- function(score, folds, seed) {
  label <- ifelse(score >= 0, 1, -1)
  assignment <- integer(length(score))
  set.seed(seed)
  counts <- table(factor(label, levels = c(-1, 1)))
  if (any(counts < folds)) {
    index <- sample.int(length(score))
    assignment[index] <- rep(seq_len(folds), length.out = length(score))
    return(assignment)
  }
  for (value in c(-1, 1)) {
    index <- sample(which(label == value))
    if (length(index) < folds) {
      .odrl_abort("Each score-sign class must contain at least `folds` rows.")
    }
    assignment[index] <- rep(seq_len(folds), length.out = length(index))
  }
  assignment
}

.odrl_log1pexp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

.odrl_hardtanh <- function(x) pmax(-1, pmin(1, x))

.odrl_empirical_criterion <- function(score, action) {
  mean(as.numeric(score) * as.numeric(action))
}

.odrl_format_seconds <- function(seconds) {
  if (!is.finite(seconds)) return("unknown")
  if (seconds < 60) return(sprintf("%.2f s", seconds))
  if (seconds < 3600) return(sprintf("%.2f min", seconds / 60))
  sprintf("%.2f h", seconds / 3600)
}
