# odrlITR

<!-- badges: start -->
[![R-CMD-check](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`odrlITR` implements **orthogonal double residual learning** for binary
individualized treatment rules.

The package separates two stages:

1. estimate `pi(x) = P(A = +1 | X = x)` and `m(x) = E(Y | X = x)` by
   cross-fitted Super Learner or built-in parametric models, or supply aligned
   nuisance predictions;
2. optimize a policy using the signed score
   `(A - e_hat(X)) * (Y - m_hat(X))`, where `e_hat = 2 * pi_hat - 1`.

## Installation

```r
# From the private GitHub repository after authenticating with GitHub:
# install.packages("pak")
# pak::pak("deckardt98/odrlITR")

# From a downloaded CRAN-style source archive:
install.packages("odrlITR_0.1.0.tar.gz", repos = NULL, type = "source")
```

Optional engines are installed only when needed:

```r
install.packages(c("SuperLearner", "glmnet", "policytree", "highs"))
```

## Quick start

### Built-in cross-fitted Super Learner

```r
library(odrlITR)

fit <- odrl(
  x = x,
  a = treatment,
  y = outcome,
  learner = "tree",
  nuisance_folds = 5,
  sl.library = c("SL.mean", "SL.glm", "SL.glmnet")
)

predict(fit, new_x)
summary(fit)
```

The outer fold count is configurable. When `nuisance_fold_id` is omitted,
`odrlITR` constructs reproducible treatment-stratified folds whose total sizes
differ by at most one row. Supply a row-aligned vector to use an externally
defined cross-fitting scheme. Every resulting training set must contain at
least two rows and, when the propensity is estimated, at least two rows from
each treatment arm:

```r
fit <- odrl(
  x, treatment, outcome,
  learner = "tree",
  nuisance_folds = 3,
  nuisance_fold_id = study_fold,
  sl.library = c("SL.mean", "SL.glm", "SL.glmnet")
)
```

Propensity and outcome libraries can be chosen separately with
`sl.library.pi` and `sl.library.m`. In a randomized trial, pass `known_pi`
(or `known_e`) to retain the known assignment mechanism while still
cross-fitting `m`.

For factor or character treatment, use `positive = "treated"` to make the
`+1` action explicit.

### Built-in parametric nuisances

Use `nuisance = "parametric"` (or `"glm"`) for cross-fitted main-effects
logistic propensity regression and Gaussian linear outcome regression. The
same configurable, treatment-stratified outer-fold mechanism is used.

```r
fit <- odrl(
  x, treatment, outcome,
  learner = "linear",
  nuisance = "parametric",
  nuisance_folds = 4
)
```

Interactions, transformations, or other prespecified basis terms can be
included directly as columns of `x`.

### User-supplied nuisances

```r
nuisance <- odrl_nuisance_user(
  m = m_oof,
  pi = pi_oof,
  fold_id = fold_id,
  source = "external cross-fitting",
  out_of_fold = TRUE
)

fit <- odrl(x, treatment, outcome, learner = "svm", loss = "logistic",
            nuisance = nuisance)
```

The nuisance argument may also be a list or a function returning a nuisance
object.

## Policy learners

| `learner` | allowed `loss` | policy class | optional package |
|---|---|---|---|
| `"tree"` | `"exact"` | shallow axis-aligned tree using `policytree` or a custom backend | `policytree` |
| `"linear"` | `"exact"` | bounded-margin affine sign rule via mixed-integer optimization | `highs` |
| `"svm"` | hinge, exponential, logistic, squared hinge, or custom | linear, Gaussian, polynomial, or custom RKHS kernel | none |
| `"relu"` | hinge, exponential, logistic, squared hinge, or custom | configurable feed-forward network or external neural backend | none |

The bounded-hinge SVM uses a normalized Gaussian kernel and an RKHS radius at
most one, which certifies scores in `[-1, 1]`. Ordinary regularized hinge fits
are also available for every supported kernel. The neural hinge learner
applies a hard-tanh output map. Exponential, logistic, squared-hinge, and custom
differentiable signed-margin losses do not have the bounded-hinge objective's
special universal-orthogonality guarantee. For signed margin \(t\),
`loss = "exponential"` uses \(\exp(-t)\). A continuously differentiable convex
linear continuation below the extreme margin \(t=-30\) prevents numerical
overflow without affecting ordinary fitted margins.

Gaussian-kernel learners form and retain an \(n \times n\) training kernel
matrix. Their time and memory costs therefore grow at least quadratically in
the training sample size, and prediction requires kernel evaluations against
all training observations. For large samples, prefer the affine, tree, or
ReLU learner, or benchmark the kernel learner on a representative subset
before committing to a full tuning grid.

Affine optimization is genuinely mixed-integer. A time-limited feasible
incumbent is not described as proven optimal: solver status, MIP gap,
constraint checks, and the objective audit are retained in
`fit$policy$diagnostics`. Time, relative-gap, absolute-gap, node, objective,
thread, logging, and additional HiGHS controls are configurable. Set
`linear_require_gap = TRUE` when all requested gap criteria must be certified.

Meeting nonzero requested gaps does not by itself prove global optimality. For
a numerical global-optimum certificate within the bounded-margin standardized
affine class, set `linear_relative_gap = 0` and leave
`linear_absolute_gap = NULL` (or set it to zero); inspect
`certified_global_optimum` in the diagnostics.

For numerical conditioning, policy optimization and tuning use the score
divided by its sample mean absolute value. This positive rescaling leaves the
unpenalized binary objective and treatment rule unchanged and makes the
package's penalty grids interpretable on a common scale. The reported training
criterion is recomputed from the original, unscaled score.

```r
control <- odrl_control(
  tree_depth = 2,
  tree_backend = "policytree",
  linear_time_limit = 300,
  linear_relative_gap = 0.01,
  linear_node_limit = 100000,
  svm_kernel = "polynomial",
  svm_polynomial_degree = c(2, 3),
  svm_hinge_mode = "regularized",
  svm_rbf_multiplier = c(0.5, 1, 2),
  svm_penalty = c(0.01, 0.1, 1),
  relu_architectures = list(integer(), 8L, c(16L, 8L)),
  relu_activation = c("relu", "tanh"),
  relu_decay = c(0.001, 0.01, 0.1),
  relu_selection = "one_se",
  seed = 2026
)
```

Custom kernels and differentiable margin losses use explicit callbacks:

```r
laplacian <- function(x, y, rate = 1) {
  distance <- sqrt(pmax(
    outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y),
    0
  ))
  exp(-rate * distance)
}

smooth_loss <- list(
  name = "softplus_margin",
  value = function(margin) {
    z <- 1 - margin
    pmax(z, 0) + log1p(exp(-abs(z)))
  },
  gradient = function(margin) -plogis(1 - margin)
)

fit <- odrl(
  x, treatment, outcome,
  learner = "svm", loss = smooth_loss,
  control = odrl_control(
    svm_kernel = list(
      name = "laplacian", fun = laplacian, args = list(rate = 0.75)
    )
  )
)
```

Custom-kernel authors are responsible for returning a symmetric
positive-semidefinite kernel; custom-loss authors are responsible for the
loss's mathematical properties, including convexity when a convex RKHS
problem is intended.

For neural policies, `relu_architectures` accepts an affine candidate
(`integer()`), one layer (`8L`), or deeper architectures such as
`c(16L, 8L)`. Built-in activations are ReLU, leaky ReLU, tanh, sigmoid, and
linear. `relu_backend` may instead provide named `fit` and `predict` callbacks
to interface with another neural-network package. Likewise, `tree_backend`
may provide callbacks for another axis-aligned tree implementation.

## Reproducibility and scope

- Outer nuisance cross-fitting and second-stage policy tuning use separate
  folds.
- Automatically generated nuisance folds are treatment-stratified and as
  nearly equal in total size as integer constraints permit.
- Treatment ties and zero scores use the paper's `+1` convention.
- Raw and numerically thresholded scores are both retained for auditing.
- The package learns a policy; causal interpretation still requires the usual
  identification conditions and appropriate nuisance estimation.
- The native neural engine supports affine, shallow, and multilayer candidates
  with configurable activations. Its hard-tanh hinge output has the paper's
  bounded-score orthogonality property, but the package does not claim the
  manuscript's sparse deep-network rate for this particular optimizer.

See `vignette("getting-started", package = "odrlITR")` for a complete workflow.

## Development status

Version 0.1.0 is an initial research-software release under private review and
being prepared for eventual CRAN submission.
