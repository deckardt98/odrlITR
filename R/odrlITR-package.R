#' Learn individualized treatment rules with ODRL
#'
#' The `odrlITR` package learns which of two treatments to recommend from
#' covariates, observed treatments, and outcomes. Larger outcomes are preferred.
#' Start with [odrl()] for a matrix or data frame, or [odrl_formula()] for a
#' formula. Use [predict.odrl_fit()] to recommend treatments for new rows.
#'
#' Orthogonal double residual learning first estimates the treatment and
#' outcome means given the covariates, then forms the training score
#' \deqn{\widehat Z_i = \{A_i-\widehat e(X_i)\}
#'       \{Y_i-\widehat m(X_i)\},}
#' where treatment is internally coded as \eqn{A\in\{-1,+1\}},
#' \eqn{e(x)=E(A\mid X=x)}, and \eqn{m(x)=E(Y\mid X=x)}.
#' Cross-fitting predicts each row from nuisance models fitted on other rows.
#' The score then guides fitting a tree, an affine rule, a kernel or finite
#' series score, or a neural network.
#'
#' The examples in `vignette("getting-started", package = "odrlITR")` use
#' `learner = "svm"` and `nuisance = "parametric"`, which need no optional
#' modeling packages. See [odrl_control()] to adjust the rule's complexity.
#'
#' @references
#' Jiaqi Tong and Fan Li (2026). "Orthogonal double residual learning for
#' optimal individualized treatment rules."
#' [arXiv:2608.24085](https://arxiv.org/abs/2608.24085).
#'
#' @keywords internal
"_PACKAGE"
