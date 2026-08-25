#' Specify a finite series score class
#'
#' Construct an explicit finite-dimensional feature map (equivalently, a
#' positive semidefinite kernel of finite rank) and fit the surrogate in the
#' primal, avoiding a dense training Gram matrix. Character shortcuts such as
#' `"legendre"`, `"fourier"`, `"bspline"`, `"haar"`, and
#' `"local_polynomial"` use the same defaults.
#'
#' Multivariate series are additive by default and therefore do not contain
#' interactions. Low-order ANOVA, total-degree Legendre, and full tensor-product
#' constructions are optional; the full tensor can grow exponentially. Bounds,
#' knots, centering, and scaling are estimated within each training fold used
#' for policy fitting, then reused for validation and prediction. Fourier terms
#' impose periodic endpoint behavior. The included wavelet is Haar, not a
#' Cohen--Daubechies--Vial boundary-corrected wavelet. The local polynomial
#' option is a partition series basis, not a local polynomial regression
#' smoother.
#'
#' @param basis Univariate basis family.
#' @param legendre_degree Positive maximum degree grid for Legendre series.
#' @param fourier_harmonics Positive maximum harmonic grid for Fourier series.
#' @param spline_df Positive grid of degrees of freedom for B-spline series.
#' @param spline_degree Nonnegative B-spline polynomial degree.
#' @param wavelet_level Positive maximum resolution grid for the Haar basis.
#' @param local_partitions Positive partition count grid for local polynomial
#'   series.
#' @param local_degree Nonnegative local polynomial degree grid.
#' @param combine Multivariate construction. `"additive"` concatenates
#'   univariate terms, `"anova"` includes products involving at most
#'   `interaction_order` variables, `"total_degree"` constructs a multivariate
#'   Legendre basis of total degree, and `"tensor"` constructs the full tensor
#'   product.
#' @param interaction_order Maximum number of variables in an ANOVA product.
#' @param variables Optional positive column indices or column names to use.
#' @param include_linear Whether to append a raw linear term for each selected
#'   coordinate when it is not already present in the Legendre basis.
#' @param domain Treatment of the coordinate domain. `"empirical"` learns
#'   bounds from each training fold; `"unit"` uses `[0,1]`.
#' @param bounds Optional fixed bounds. Supply a numeric vector of length two,
#'   a list with `lower` and `upper`, or a matrix with two columns.
#' @param boundary_quantiles Quantiles from each training fold used for
#'   empirical bounds.
#' @param extrapolation Whether prediction values outside fitted bounds are
#'   clamped or rejected.
#' @param normalize Whether to divide centered series columns by their
#'   root-mean-square values from each training fold.
#' @param max_features Maximum feature count before centering and rank removal.
#' @param max_feature_elements Maximum product of rows and generated features.
#'
#' @return An object of class `odrl_series_kernel`.
#' @export
#'
#' @examples
#' series <- odrl_series_kernel(
#'   basis = "legendre", legendre_degree = c(1, 2, 3),
#'   combine = "anova", interaction_order = 2
#' )
#' control <- odrl_control(
#'   svm_kernel = series, svm_penalty = c(0.01, 0.1, 1)
#' )
odrl_series_kernel <- function(
    basis = c("legendre", "fourier", "bspline", "haar",
              "local_polynomial"),
    legendre_degree = c(1L, 2L, 3L),
    fourier_harmonics = c(1L, 2L, 3L),
    spline_df = c(4L, 6L, 8L),
    spline_degree = 3L,
    wavelet_level = c(1L, 2L, 3L),
    local_partitions = c(2L, 4L),
    local_degree = c(0L, 1L),
    combine = c("additive", "anova", "total_degree", "tensor"),
    interaction_order = 1L,
    variables = NULL,
    include_linear = FALSE,
    domain = c("empirical", "unit"),
    bounds = NULL,
    boundary_quantiles = c(0.01, 0.99),
    extrapolation = c("clamp", "error"),
    normalize = TRUE,
    max_features = 1000L,
    max_feature_elements = 2.5e7) {
  if (!is.character(basis) || !length(basis) || anyNA(basis)) {
    .odrl_abort("`basis` must contain a supported series basis name.")
  }
  basis <- .odrl_series_basis_name(basis[[1L]])
  combine <- match.arg(combine)
  domain <- match.arg(domain)
  extrapolation <- match.arg(extrapolation)
  legendre_degree <- .odrl_series_integer_grid(
    legendre_degree, "legendre_degree", 1L
  )
  fourier_harmonics <- .odrl_series_integer_grid(
    fourier_harmonics, "fourier_harmonics", 1L
  )
  spline_df <- .odrl_series_integer_grid(spline_df, "spline_df", 1L)
  .odrl_check_scalar(spline_degree, "spline_degree", 0, Inf, integer = TRUE)
  if (any(spline_df < as.integer(spline_degree) + 1L)) {
    .odrl_abort(
      "Every `spline_df` value must be at least `spline_degree + 1`."
    )
  }
  wavelet_level <- .odrl_series_integer_grid(
    wavelet_level, "wavelet_level", 1L
  )
  local_partitions <- .odrl_series_integer_grid(
    local_partitions, "local_partitions", 2L
  )
  local_degree <- .odrl_series_integer_grid(
    local_degree, "local_degree", 0L
  )
  .odrl_check_scalar(
    interaction_order, "interaction_order", 1, Inf, integer = TRUE
  )
  if (combine == "total_degree" && basis != "legendre") {
    .odrl_abort("`combine = \"total_degree\"` requires a Legendre basis.")
  }
  if (combine == "total_degree" && isTRUE(include_linear)) {
    .odrl_abort(
      "`include_linear` is redundant for a Legendre basis of total degree."
    )
  }
  if (!is.null(variables)) {
    valid_numeric <- is.numeric(variables) && length(variables) &&
      all(is.finite(variables)) && all(variables >= 1) &&
      all(variables == as.integer(variables))
    valid_character <- is.character(variables) && length(variables) &&
      !anyNA(variables) && all(nzchar(variables))
    if (!valid_numeric && !valid_character) {
      .odrl_abort(
        "`variables` must be NULL, positive column indices, or column names."
      )
    }
    if (anyDuplicated(variables)) {
      .odrl_abort("`variables` must not contain duplicates.")
    }
  }
  .odrl_series_flag(include_linear, "include_linear")
  .odrl_series_flag(normalize, "normalize")
  if (!is.numeric(boundary_quantiles) || length(boundary_quantiles) != 2L ||
      any(!is.finite(boundary_quantiles)) ||
      boundary_quantiles[[1L]] < 0 || boundary_quantiles[[2L]] > 1 ||
      boundary_quantiles[[1L]] >= boundary_quantiles[[2L]]) {
    .odrl_abort(
      "`boundary_quantiles` must be two increasing numbers in [0,1]."
    )
  }
  .odrl_series_validate_bounds(bounds)
  .odrl_check_scalar(max_features, "max_features", 1, Inf, integer = TRUE)
  .odrl_check_scalar(
    max_feature_elements, "max_feature_elements", 1, Inf
  )
  structure(list(
    name = basis,
    basis = basis,
    legendre_degree = legendre_degree,
    fourier_harmonics = fourier_harmonics,
    spline_df = spline_df,
    spline_degree = as.integer(spline_degree),
    wavelet_level = wavelet_level,
    local_partitions = local_partitions,
    local_degree = local_degree,
    combine = combine,
    interaction_order = as.integer(interaction_order),
    variables = variables,
    include_linear = include_linear,
    domain = domain,
    bounds = bounds,
    boundary_quantiles = as.numeric(boundary_quantiles),
    extrapolation = extrapolation,
    normalize = normalize,
    max_features = as.integer(max_features),
    max_feature_elements = as.numeric(max_feature_elements),
    builtin = TRUE,
    finite_features = TRUE
  ), class = c("odrl_series_kernel", "list"))
}

.odrl_series_basis_name <- function(basis) {
  if (!is.character(basis) || length(basis) != 1L || is.na(basis)) {
    .odrl_abort("A series basis name must be one nonmissing string.")
  }
  value <- tolower(gsub("-", "_", basis, fixed = TRUE))
  value <- switch(value,
    b_spline = "bspline",
    spline = "bspline",
    wavelet = "haar",
    wavelet_haar = "haar",
    local_poly = "local_polynomial",
    partition = "local_polynomial",
    value
  )
  allowed <- c("legendre", "fourier", "bspline", "haar",
               "local_polynomial")
  if (!value %in% allowed) {
    .odrl_abort(
      "Unknown series basis `", basis, "`. Use `\"legendre\"`, ",
      "`\"fourier\"`, `\"bspline\"`, `\"haar\"`, or ",
      "`\"local_polynomial\"`."
    )
  }
  value
}

.odrl_series_integer_grid <- function(value, name, minimum) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
      any(value < minimum) || any(value != as.integer(value))) {
    .odrl_abort("`", name, "` must contain integers >= ", minimum, ".")
  }
  sort(unique(as.integer(value)))
}

.odrl_series_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    .odrl_abort("`", name, "` must be TRUE or FALSE.")
  }
  invisible(value)
}

.odrl_series_validate_bounds <- function(bounds) {
  if (is.null(bounds)) return(invisible(NULL))
  if (is.numeric(bounds) && is.null(dim(bounds))) {
    if (length(bounds) != 2L || any(!is.finite(bounds)) ||
        bounds[[1L]] >= bounds[[2L]]) {
      .odrl_abort("Numeric `bounds` must be a finite increasing pair.")
    }
    return(invisible(NULL))
  }
  if (is.list(bounds)) {
    if (!all(c("lower", "upper") %in% names(bounds)) ||
        !is.numeric(bounds$lower) || !is.numeric(bounds$upper) ||
        any(!is.finite(bounds$lower)) || any(!is.finite(bounds$upper))) {
      .odrl_abort("List `bounds` must contain finite `lower` and `upper`.")
    }
    return(invisible(NULL))
  }
  if (is.matrix(bounds) && is.numeric(bounds) && ncol(bounds) == 2L &&
      all(is.finite(bounds))) {
    return(invisible(NULL))
  }
  .odrl_abort(
    "`bounds` must be NULL, an increasing numeric pair, a `lower`/`upper` ",
    "list, or a finite two-column matrix."
  )
}

.odrl_is_series_kernel <- function(kernel) {
  if (inherits(kernel, "odrl_series_kernel")) return(TRUE)
  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel)) {
    return(FALSE)
  }
  value <- tolower(gsub("-", "_", kernel, fixed = TRUE))
  value %in% c(
    "legendre", "fourier", "bspline", "b_spline", "spline", "haar",
    "wavelet", "wavelet_haar", "local_polynomial", "local_poly",
    "partition"
  )
}

.odrl_resolve_series_kernel <- function(kernel) {
  if (inherits(kernel, "odrl_series_kernel")) return(kernel)
  if (.odrl_is_series_kernel(kernel)) return(odrl_series_kernel(kernel))
  .odrl_abort("The supplied object is not a series kernel specification.")
}

.odrl_series_grid <- function(spec) {
  spec <- .odrl_resolve_series_kernel(spec)
  switch(spec$basis,
    legendre = data.frame(degree = spec$legendre_degree),
    fourier = data.frame(harmonics = spec$fourier_harmonics),
    bspline = data.frame(df = spec$spline_df),
    haar = data.frame(level = spec$wavelet_level),
    local_polynomial = expand.grid(
      partitions = spec$local_partitions,
      local_degree = spec$local_degree,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

.odrl_series_candidate <- function(spec, candidate = list()) {
  spec <- .odrl_resolve_series_kernel(spec)
  first <- .odrl_series_grid(spec)[1L, , drop = FALSE]
  value <- as.list(first)
  for (name in intersect(names(candidate), names(value))) {
    if (length(candidate[[name]]) == 1L && !is.na(candidate[[name]])) {
      value[[name]] <- candidate[[name]]
    }
  }
  value
}

.odrl_series_matrix <- function(x, name = "x") {
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x))) {
    .odrl_abort("`", name, "` must be a nonempty finite numeric matrix.")
  }
  x
}

.odrl_series_variables <- function(x, variables) {
  if (is.null(variables)) return(seq_len(ncol(x)))
  if (is.character(variables)) {
    if (is.null(colnames(x))) {
      .odrl_abort("Character `variables` require column names in `x`.")
    }
    index <- match(variables, colnames(x))
    if (anyNA(index)) {
      .odrl_abort("Unknown series variable(s): ",
                  paste(variables[is.na(index)], collapse = ", "), ".")
    }
    return(as.integer(index))
  }
  index <- as.integer(variables)
  if (any(index > ncol(x))) {
    .odrl_abort("A series variable index exceeds `ncol(x)`.")
  }
  index
}

.odrl_series_expand_bound <- function(value, p, name) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
      !length(value) %in% c(1L, p)) {
    .odrl_abort("`bounds$", name, "` must have length one or ", p, ".")
  }
  rep(as.numeric(value), length.out = p)
}

.odrl_series_bounds <- function(x, spec) {
  p <- ncol(x)
  if (!is.null(spec$bounds)) {
    bounds <- spec$bounds
    if (is.numeric(bounds) && is.null(dim(bounds))) {
      lower <- rep(bounds[[1L]], p)
      upper <- rep(bounds[[2L]], p)
    } else if (is.list(bounds)) {
      lower <- .odrl_series_expand_bound(bounds$lower, p, "lower")
      upper <- .odrl_series_expand_bound(bounds$upper, p, "upper")
    } else {
      if (nrow(bounds) == 1L) bounds <- bounds[rep(1L, p), , drop = FALSE]
      if (nrow(bounds) != p) {
        .odrl_abort("Matrix `bounds` must have one row per selected variable.")
      }
      lower <- bounds[, 1L]
      upper <- bounds[, 2L]
    }
  } else if (identical(spec$domain, "unit")) {
    lower <- rep(0, p)
    upper <- rep(1, p)
  } else {
    probability <- spec$boundary_quantiles
    lower <- apply(x, 2L, stats::quantile, probs = probability[[1L]],
                   names = FALSE, type = 8)
    upper <- apply(x, 2L, stats::quantile, probs = probability[[2L]],
                   names = FALSE, type = 8)
    constant <- upper - lower <= sqrt(.Machine$double.eps) *
      pmax(1, abs(lower), abs(upper))
    if (any(constant)) {
      midpoint <- (lower[constant] + upper[constant]) / 2
      lower[constant] <- midpoint - 0.5
      upper[constant] <- midpoint + 0.5
    }
  }
  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
      any(lower >= upper)) {
    .odrl_abort("Every fitted series bound must be finite and increasing.")
  }
  list(lower = as.numeric(lower), upper = as.numeric(upper))
}

.odrl_series_apply_bounds <- function(x, bounds,
                                      extrapolation = "clamp") {
  tolerance <- 1e-12 * pmax(1, abs(bounds$lower), abs(bounds$upper))
  below <- sweep(x, 2L, bounds$lower - tolerance, "<")
  above <- sweep(x, 2L, bounds$upper + tolerance, ">")
  if (identical(extrapolation, "error") && any(below | above)) {
    .odrl_abort("New data lie outside the fitted series bounds.")
  }
  u <- sweep(x, 2L, bounds$lower, "-")
  u <- sweep(u, 2L, bounds$upper - bounds$lower, "/")
  list(
    u = pmin(pmax(u, 0), 1),
    clamped_below = colSums(below),
    clamped_above = colSums(above)
  )
}

.odrl_legendre_matrix <- function(u, degree) {
  t <- 2 * u - 1
  result <- matrix(0, nrow = length(u), ncol = degree)
  previous <- rep(1, length(u))
  current <- t
  result[, 1L] <- sqrt(3) * current
  if (degree >= 2L) {
    for (order in 2:degree) {
      next_value <- ((2 * order - 1) * t * current -
        (order - 1) * previous) / order
      result[, order] <- sqrt(2 * order + 1) * next_value
      previous <- current
      current <- next_value
    }
  }
  colnames(result) <- paste0("P", seq_len(degree))
  result
}

.odrl_fourier_matrix <- function(u, harmonics) {
  result <- matrix(0, nrow = length(u), ncol = 2L * harmonics)
  names <- character(2L * harmonics)
  for (harmonic in seq_len(harmonics)) {
    cosine <- 2L * harmonic - 1L
    sine <- 2L * harmonic
    result[, cosine] <- sqrt(2) * cos(2 * pi * harmonic * u)
    result[, sine] <- sqrt(2) * sin(2 * pi * harmonic * u)
    names[c(cosine, sine)] <- c(
      paste0("cos", harmonic), paste0("sin", harmonic)
    )
  }
  colnames(result) <- names
  result
}

.odrl_haar_matrix <- function(u, level) {
  count <- 2^level - 1L
  result <- matrix(0, nrow = length(u), ncol = count)
  names <- character(count)
  adjusted <- pmin(u, 1 - .Machine$double.eps)
  cursor <- 0L
  for (resolution in 0:(level - 1L)) {
    cells <- 2^resolution
    position <- adjusted * cells
    cell <- floor(position)
    sign <- ifelse(position - cell < 0.5, 1, -1)
    for (location in 0:(cells - 1L)) {
      cursor <- cursor + 1L
      result[, cursor] <- sqrt(cells) * sign * (cell == location)
      names[[cursor]] <- paste0("psi", resolution, "_", location)
    }
  }
  colnames(result) <- names
  result
}

.odrl_local_polynomial_matrix <- function(u, partitions, degree) {
  adjusted <- pmin(u, 1 - .Machine$double.eps)
  bin <- floor(adjusted * partitions)
  local <- adjusted * partitions - (bin + 0.5)
  count <- partitions * (degree + 1L) - 1L
  result <- matrix(0, nrow = length(u), ncol = count)
  names <- character(count)
  cursor <- 0L
  for (order in 0:degree) {
    for (cell in 0:(partitions - 1L)) {
      if (order == 0L && cell == partitions - 1L) next
      cursor <- cursor + 1L
      result[, cursor] <- (bin == cell) * local^order
      names[[cursor]] <- paste0("bin", cell + 1L, "_d", order)
    }
  }
  colnames(result) <- names
  result
}

.odrl_series_univariate_fit <- function(u, spec, candidate) {
  basis <- spec$basis
  state <- list(basis = basis, include_linear = FALSE)
  result <- switch(basis,
    legendre = {
      state$degree <- as.integer(candidate$degree)
      .odrl_legendre_matrix(u, state$degree)
    },
    fourier = {
      state$harmonics <- as.integer(candidate$harmonics)
      .odrl_fourier_matrix(u, state$harmonics)
    },
    bspline = {
      state$df <- as.integer(candidate$df)
      state$degree <- spec$spline_degree
      full <- splines::bs(
        u, df = state$df, degree = state$degree, intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
      state$knots <- attr(full, "knots")
      full <- full[, -ncol(full), drop = FALSE]
      colnames(full) <- paste0("B", seq_len(ncol(full)))
      full
    },
    haar = {
      state$level <- as.integer(candidate$level)
      .odrl_haar_matrix(u, state$level)
    },
    local_polynomial = {
      state$partitions <- as.integer(candidate$partitions)
      state$degree <- as.integer(candidate$local_degree)
      .odrl_local_polynomial_matrix(u, state$partitions, state$degree)
    }
  )
  if (isTRUE(spec$include_linear) && basis != "legendre") {
    result <- cbind(result, linear = u - 0.5)
    state$include_linear <- TRUE
  }
  list(features = result, state = state)
}

.odrl_series_univariate_apply <- function(u, state) {
  result <- switch(state$basis,
    legendre = .odrl_legendre_matrix(u, state$degree),
    fourier = .odrl_fourier_matrix(u, state$harmonics),
    bspline = {
      full <- splines::bs(
        u, knots = state$knots, degree = state$degree, intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
      full <- full[, -ncol(full), drop = FALSE]
      colnames(full) <- paste0("B", seq_len(ncol(full)))
      full
    },
    haar = .odrl_haar_matrix(u, state$level),
    local_polynomial = .odrl_local_polynomial_matrix(
      u, state$partitions, state$degree
    )
  )
  if (isTRUE(state$include_linear)) result <- cbind(result, linear = u - 0.5)
  result
}

.odrl_series_feature_count <- function(spec, candidate, p) {
  univariate <- switch(spec$basis,
    legendre = as.numeric(candidate$degree),
    fourier = 2 * as.numeric(candidate$harmonics),
    bspline = as.numeric(candidate$df) - 1,
    haar = 2^as.numeric(candidate$level) - 1,
    local_polynomial = as.numeric(candidate$partitions) *
      (as.numeric(candidate$local_degree) + 1) - 1
  )
  if (isTRUE(spec$include_linear) && spec$basis != "legendre") {
    univariate <- univariate + 1
  }
  count <- switch(spec$combine,
    additive = p * univariate,
    anova = {
      order <- min(spec$interaction_order, p)
      sum(vapply(seq_len(order), function(size) {
        choose(p, size) * univariate^size
      }, numeric(1)))
    },
    tensor = (univariate + 1)^p - 1,
    total_degree = choose(p + as.numeric(candidate$degree),
                          as.numeric(candidate$degree)) - 1
  )
  as.numeric(count)
}

.odrl_series_guard_features <- function(spec, candidate, p, n) {
  count <- .odrl_series_feature_count(spec, candidate, p)
  elements <- n * count
  if (!is.finite(count) || count > spec$max_features) {
    .odrl_abort(
      "The requested series contains ", format(count, scientific = FALSE),
      " features, exceeding `max_features = ", spec$max_features, "`."
    )
  }
  if (!is.finite(elements) || elements > spec$max_feature_elements) {
    .odrl_abort(
      "The requested series design contains ",
      format(elements, scientific = FALSE),
      " elements, exceeding `max_feature_elements = ",
      format(spec$max_feature_elements, scientific = FALSE), "`."
    )
  }
  as.integer(round(count))
}

.odrl_series_row_tensor <- function(blocks) {
  result <- blocks[[1L]]
  if (length(blocks) == 1L) return(result)
  for (block in blocks[-1L]) {
    left <- result
    result <- matrix(0, nrow = nrow(left),
                     ncol = ncol(left) * ncol(block))
    names <- character(ncol(result))
    cursor <- 0L
    for (i in seq_len(ncol(left))) {
      for (j in seq_len(ncol(block))) {
        cursor <- cursor + 1L
        result[, cursor] <- left[, i] * block[, j]
        names[[cursor]] <- paste(colnames(left)[[i]],
                                 colnames(block)[[j]], sep = ":")
      }
    }
    colnames(result) <- names
  }
  result
}

.odrl_total_degree_indices <- function(p, degree) {
  output <- list()
  recurse <- function(position, remaining, current) {
    if (position == p) {
      current[[position]] <- remaining
      output[[length(output) + 1L]] <<- current
      return(invisible(NULL))
    }
    for (value in 0:remaining) {
      current[[position]] <- value
      recurse(position + 1L, remaining - value, current)
    }
    invisible(NULL)
  }
  for (total in seq_len(degree)) recurse(1L, total, integer(p))
  output
}

.odrl_series_combine <- function(blocks, variable_names, spec, candidate) {
  for (j in seq_along(blocks)) {
    colnames(blocks[[j]]) <- paste0(
      variable_names[[j]], "_", colnames(blocks[[j]])
    )
  }
  p <- length(blocks)
  if (spec$combine == "additive") return(do.call(cbind, blocks))
  if (spec$combine == "total_degree") {
    indices <- .odrl_total_degree_indices(p, as.integer(candidate$degree))
    result <- matrix(1, nrow = nrow(blocks[[1L]]), ncol = length(indices))
    names <- character(length(indices))
    for (column in seq_along(indices)) {
      alpha <- indices[[column]]
      pieces <- character()
      for (j in which(alpha > 0L)) {
        result[, column] <- result[, column] * blocks[[j]][, alpha[[j]]]
        pieces <- c(pieces, colnames(blocks[[j]])[[alpha[[j]]]])
      }
      names[[column]] <- paste(pieces, collapse = ":")
    }
    colnames(result) <- names
    return(result)
  }
  maximum_order <- if (spec$combine == "tensor") p else
    min(spec$interaction_order, p)
  output <- vector("list", 0L)
  cursor <- 0L
  for (order in seq_len(maximum_order)) {
    subsets <- utils::combn(seq_len(p), order, simplify = FALSE)
    for (subset in subsets) {
      cursor <- cursor + 1L
      output[[cursor]] <- .odrl_series_row_tensor(blocks[subset])
    }
  }
  do.call(cbind, output)
}

.odrl_fit_series_map <- function(x, spec, candidate = list()) {
  x <- .odrl_series_matrix(x)
  spec <- .odrl_resolve_series_kernel(spec)
  candidate <- .odrl_series_candidate(spec, candidate)
  variables <- .odrl_series_variables(x, spec$variables)
  selected <- x[, variables, drop = FALSE]
  variable_names <- colnames(x)[variables]
  if (is.null(variable_names)) variable_names <- paste0("x", variables)
  .odrl_series_guard_features(
    spec, candidate, length(variables), nrow(selected)
  )
  bounds <- .odrl_series_bounds(selected, spec)
  domain <- .odrl_series_apply_bounds(
    selected, bounds, spec$extrapolation
  )
  trained <- lapply(seq_len(ncol(domain$u)), function(j) {
    .odrl_series_univariate_fit(domain$u[, j], spec, candidate)
  })
  blocks <- lapply(trained, `[[`, "features")
  states <- lapply(trained, `[[`, "state")
  raw <- .odrl_series_combine(
    blocks, variable_names, spec, candidate
  )
  if (!is.matrix(raw) || !ncol(raw) || any(!is.finite(raw))) {
    .odrl_abort("The fitted series basis produced no finite features.")
  }
  center <- colMeans(raw)
  centered <- sweep(raw, 2L, center, "-")
  rms <- sqrt(colMeans(centered^2))
  keep <- is.finite(rms) & rms > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    .odrl_abort("Every fitted series feature is constant on the training data.")
  }
  centered <- centered[, keep, drop = FALSE]
  divisor <- if (isTRUE(spec$normalize)) rms[keep] else rep(1, sum(keep))
  feature_scale <- sqrt(sum(keep))
  features <- sweep(centered, 2L, divisor, "/") / feature_scale
  map <- structure(list(
    spec = spec,
    candidate = candidate,
    variables = variables,
    input_columns = ncol(x),
    input_names = colnames(x),
    variable_names = variable_names,
    bounds = bounds,
    univariate_state = states,
    center = center,
    rms = rms,
    keep = keep,
    divisor = divisor,
    feature_scale = feature_scale,
    term_names = colnames(raw)[keep],
    raw_feature_count = ncol(raw),
    feature_count = sum(keep),
    training_clamped_below = domain$clamped_below,
    training_clamped_above = domain$clamped_above
  ), class = c("odrl_series_map", "list"))
  colnames(features) <- map$term_names
  list(features = features, map = map)
}

.odrl_apply_series_map <- function(map, newx, diagnostics = FALSE) {
  if (!inherits(map, "odrl_series_map")) {
    .odrl_abort("`map` must be a fitted `odrl_series_map` object.")
  }
  newx <- .odrl_series_matrix(newx, "newx")
  if (ncol(newx) != map$input_columns) {
    .odrl_abort("`newx` has a different number of columns from the training data.")
  }
  selected <- newx[, map$variables, drop = FALSE]
  domain <- .odrl_series_apply_bounds(
    selected, map$bounds, map$spec$extrapolation
  )
  blocks <- lapply(seq_along(map$univariate_state), function(j) {
    .odrl_series_univariate_apply(
      domain$u[, j], map$univariate_state[[j]]
    )
  })
  raw <- .odrl_series_combine(
    blocks, map$variable_names, map$spec, map$candidate
  )
  raw <- raw[, map$keep, drop = FALSE]
  features <- sweep(raw, 2L, map$center[map$keep], "-")
  features <- sweep(features, 2L, map$divisor, "/") / map$feature_scale
  colnames(features) <- map$term_names
  if (!isTRUE(diagnostics)) return(features)
  list(
    features = features,
    clamped_below = domain$clamped_below,
    clamped_above = domain$clamped_above
  )
}

.odrl_fit_series_surrogate <- function(
    x, score, spec, candidate = list(), lambda = 0.1, maxit = 300L,
    seed = 1L, loss = "logistic") {
  x <- .odrl_series_matrix(x)
  if (!is.numeric(score) || length(score) != nrow(x) ||
      any(!is.finite(score))) {
    .odrl_abort("`score` must contain one finite number per row of `x`.")
  }
  .odrl_check_scalar(lambda, "lambda", 0, Inf, open_lower = TRUE)
  .odrl_check_scalar(maxit, "maxit", 1, Inf, integer = TRUE)
  .odrl_check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  loss_spec <- .odrl_resolve_surrogate_loss(loss)
  mapped <- .odrl_fit_series_map(x, spec, candidate)
  features <- mapped$features
  label <- ifelse(score >= 0, 1, -1)
  score_scale <- mean(abs(score))
  weight <- if (is.finite(score_scale) && score_scale > 0) {
    abs(score) / score_scale
  } else {
    rep(1, length(score))
  }
  n <- length(score)
  objective <- function(theta) {
    raw <- theta[[1L]] + drop(features %*% theta[-1L])
    margin <- label * raw
    mean(weight * loss_spec$value(margin)) +
      0.5 * lambda * sum(theta^2)
  }
  gradient <- function(theta) {
    raw <- theta[[1L]] + drop(features %*% theta[-1L])
    margin <- label * raw
    derivative <- weight * loss_spec$gradient(margin) * label / n
    c(sum(derivative), drop(crossprod(features, derivative))) + lambda * theta
  }
  initial_intercept <- if (isTRUE(loss_spec$builtin) &&
      identical(loss_spec$name, "exponential")) {
    positive <- max(sum(weight[label == 1]), .Machine$double.eps)
    negative <- max(sum(weight[label == -1]), .Machine$double.eps)
    0.5 * log(positive / negative)
  } else {
    stats::qlogis(pmin(pmax(mean(label == 1), 0.01), 0.99))
  }
  set.seed(seed)
  fit <- stats::optim(
    par = c(initial_intercept, rep(0, ncol(features))),
    fn = objective, gr = gradient, method = "L-BFGS-B",
    control = list(maxit = as.integer(maxit), factr = 1e8)
  )
  attempts <- 1L
  if (fit$convergence != 0L) {
    retry <- stats::optim(
      par = fit$par, fn = objective, gr = gradient, method = "L-BFGS-B",
      control = list(
        maxit = max(3L * as.integer(maxit), as.integer(maxit) + 100L),
        factr = 1e8
      )
    )
    attempts <- 2L
    if (retry$convergence == 0L || retry$value < fit$value) fit <- retry
  }
  fitted <- fit$par[[1L]] + drop(features %*% fit$par[-1L])
  list(
    intercept = fit$par[[1L]],
    beta = fit$par[-1L],
    series_map = mapped$map,
    kernel_spec = .odrl_resolve_series_kernel(spec),
    candidate = mapped$map$candidate,
    lambda = lambda,
    loss = loss_spec$name,
    loss_spec = loss_spec,
    fitted = fitted,
    convergence = as.integer(fit$convergence),
    message = fit$message,
    attempts = attempts,
    objective = fit$value,
    rkhs_norm = sqrt(sum(fit$par^2)),
    globally_bounded = FALSE,
    bounded_output = FALSE,
    hinge_mode = if (isTRUE(loss_spec$builtin) &&
      identical(loss_spec$name, "hinge")) "regularized" else NULL
  )
}

.odrl_predict_series_unclipped <- function(fit, newx) {
  if (is.null(fit$series_map) || is.null(fit$beta)) {
    .odrl_abort("The supplied object is not a fitted series surrogate.")
  }
  features <- .odrl_apply_series_map(fit$series_map, newx)
  fit$intercept + drop(features %*% fit$beta)
}
