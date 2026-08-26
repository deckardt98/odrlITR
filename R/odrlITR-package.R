#' Orthogonal Double Residual Learning
#'
#' The `odrlITR` package learns binary individualized treatment rules from the
#' cross-fitted signed score
#' \deqn{\widehat Z_i = \{A_i-\widehat e(X_i)\}
#'       \{Y_i-\widehat m(X_i)\},}
#' where treatment is internally coded as \eqn{A\in\{-1,+1\}},
#' \eqn{e(x)=E(A\mid X=x)}, and \eqn{m(x)=E(Y\mid X=x)}.
#'
#' @references
#' Tong, J., and Li, F. (2026). "Orthogonal double residual learning for
#' optimal individualized treatment rules."
#' [arXiv:2608.24085](https://arxiv.org/abs/2608.24085).
#'
#' @keywords internal
"_PACKAGE"
