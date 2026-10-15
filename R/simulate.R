#' Simulate data with a known treatment rule
#'
#' Generate five normally distributed covariates, binary treatment, and a
#' numeric outcome for examples and tests. Treatment follows a logistic model,
#' the marginal outcome mean is linear, and the best treatment is determined
#' by either a linear boundary or a tree with depth two.
#'
#' @param n Number of observations, at least 50.
#' @param boundary Shape of the true treatment rule, either `"linear"` or
#'   `"tree"`.
#' @param seed Random seed.
#'
#' @return A list with one row or entry per observation.
#'   \describe{
#'     \item{`x`}{A numeric matrix with five covariates, named `X1` to `X5`.}
#'     \item{`a`}{Observed treatment, coded as `-1` or `+1`.}
#'     \item{`y`}{Observed outcome, with larger values preferred.}
#'     \item{`pi`}{True probability of receiving treatment `+1`.}
#'     \item{`e`}{True conditional treatment mean, `2 * pi - 1`.}
#'     \item{`m`}{True marginal outcome mean given the covariates.}
#'     \item{`tau`}{True conditional mean outcome difference between
#'       treatments `+1` and `-1`.}
#'     \item{`optimal`}{The true optimal treatment, `+1` when `tau >= 0`
#'       and `-1` otherwise.}
#'   }
#' @export
#'
#' @examples
#' dat <- odrl_simulate(n = 100, seed = 1)
#' head(dat$x)
#' table(dat$a)
#' head(dat$optimal)
odrl_simulate <- function(n = 500L, boundary = c("linear", "tree"), seed = 1L) {
  boundary <- match.arg(boundary)
  .odrl_check_scalar(n, "n", 50, Inf, integer = TRUE)
  set.seed(seed)
  x <- matrix(stats::rnorm(n * 5L), nrow = n, ncol = 5L)
  colnames(x) <- paste0("X", seq_len(ncol(x)))
  # Treatment depends on the first two covariates.
  pi <- stats::plogis(0.4 * x[, 1L] - 0.3 * x[, 2L])
  a01 <- stats::rbinom(n, 1L, pi)
  a <- 2 * a01 - 1
  # The pooled outcome mean and treatment difference play separate roles.
  m <- 1 + x[, 3L] - 0.5 * x[, 4L]
  if (boundary == "linear") {
    tau <- x[, 1L] + x[, 2L]
  } else {
    tau <- ifelse(x[, 1L] < 0, sign(x[, 2L]), sign(x[, 3L]))
  }
  # Center treatment by its conditional mean so E(Y | X) remains m.
  y <- m + (a - (2 * pi - 1)) * tau / 2 + stats::rnorm(n)
  list(
    x = x, a = a, y = y, pi = pi, e = 2 * pi - 1, m = m,
    tau = tau, optimal = ifelse(tau >= 0, 1, -1)
  )
}
