#' Simulate a simple ODRL example
#'
#' Generates a binary-treatment example with a linear propensity, linear
#' marginal outcome regression, and either a linear or depth-2 treatment rule.
#' It is intended for examples and smoke tests, not as a scientific benchmark.
#'
#' @param n Sample size.
#' @param boundary Either `"linear"` or `"tree"`.
#' @param seed Random seed.
#'
#' @return A list containing `x`, `a`, `y`, and the true nuisances and rule.
#' @export
odrl_simulate <- function(n = 500L, boundary = c("linear", "tree"), seed = 1L) {
  boundary <- match.arg(boundary)
  .odrl_check_scalar(n, "n", 50, Inf, integer = TRUE)
  set.seed(seed)
  x <- matrix(stats::rnorm(n * 5L), nrow = n, ncol = 5L)
  colnames(x) <- paste0("X", seq_len(ncol(x)))
  pi <- stats::plogis(0.4 * x[, 1L] - 0.3 * x[, 2L])
  a01 <- stats::rbinom(n, 1L, pi)
  a <- 2 * a01 - 1
  m <- 1 + x[, 3L] - 0.5 * x[, 4L]
  if (boundary == "linear") {
    tau <- x[, 1L] + x[, 2L]
  } else {
    tau <- ifelse(x[, 1L] < 0, sign(x[, 2L]), sign(x[, 3L]))
  }
  y <- m + (a - (2 * pi - 1)) * tau / 2 + stats::rnorm(n)
  list(
    x = x, a = a, y = y, pi = pi, e = 2 * pi - 1, m = m,
    tau = tau, optimal = ifelse(tau >= 0, 1, -1)
  )
}
