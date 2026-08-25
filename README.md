# odrlITR

<!-- badges: start -->
[![R-CMD-check](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`odrlITR` implements **orthogonal double residual learning** for estimating optimal
individualized treatment rules.

The package implements a two-stage algorithm:

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
install.packages(c("SuperLearner", "glmnet", "policytree", "highs", "nnet"))
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
logistic propensity regression and linear outcome regression. The
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

The interface follows a quick-start-plus-override design: each policy class has
small built-in defaults and presets, while kernels, losses, fold assignments,
tree engines, and neural engines can be supplied by the user.

| `learner` | allowed `loss` | policy class | optional package |
|---|---|---|---|
| `"tree"` | `"exact"` | shallow axis-aligned tree using `policytree` or a custom backend | `policytree` |
| `"linear"` | `"exact"` | bounded-margin affine sign rule via mixed-integer optimization | `highs` |
| `"svm"` | hinge, exponential, logistic, squared hinge, or custom | linear, Gaussian, polynomial, finite-series, or custom RKHS score class | none |
| `"relu"` | hinge, exponential, logistic, squared hinge, or custom | configurable feed-forward network or external neural backend | `nnet` for its optional logistic-only backend |

The bounded-hinge SVM uses a normalized Gaussian kernel and an RKHS radius at
most one, which certifies scores in `[-1, 1]`. Ordinary regularized hinge fits
are also available for every supported kernel. The neural hinge learner
applies a hard-tanh output map. Exponential, logistic, squared-hinge, and custom
differentiable signed-margin losses do not have the bounded-hinge objective's
special universal-orthogonality guarantee. For signed margin \(t\),
`loss = "exponential"` uses \(\exp(-t)\). A continuously differentiable convex
linear continuation below the extreme margin \(t=-30\) prevents numerical
overflow without affecting ordinary fitted margins.

Gaussian and general custom-kernel learners form and retain an \(n \times n\) training kernel
matrix. Their time and memory costs therefore grow at least quadratically in
the training sample size, and prediction requires kernel evaluations against
all training observations. Finite-series learners instead use an explicit
feature map and a primal fit, avoiding this dense kernel matrix. For large
samples, prefer an affine, tree, finite-series, or neural learner, or benchmark
a dense kernel learner before committing to a full tuning grid.

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

For an exact tree, `tree_depth` is the **maximum** depth: zero gives a constant
rule, one a stump, and two permits a root split followed by child splits. The
fit need not use every allowed split. Search cost can rise rapidly with depth;
`tree_split_step > 1` reduces computation by thinning candidate split points
and therefore changes the searched policy class.

### Finite-series score classes

The package provides leakage-safe, training-fold-fitted feature maps for common
series approximations:

| shortcut | basis | useful controls |
|---|---|---|
| `"legendre"` | Legendre polynomial series | `legendre_degree` |
| `"fourier"` | sine/cosine harmonics | `fourier_harmonics` |
| `"bspline"` | B-spline series | `spline_df`, `spline_degree` |
| `"haar"` | Haar wavelets | `wavelet_level` |
| `"local_polynomial"` | partition indicators and centered powers | `local_partitions`, `local_degree` |

```r
series <- odrl_series_kernel(
  basis = "legendre",
  legendre_degree = c(1, 2, 3),
  combine = "anova",
  interaction_order = 2,
  max_features = 1000
)

fit <- odrl(
  x, treatment, outcome,
  learner = "svm", loss = "logistic",
  control = odrl_control(
    svm_kernel = series,
    svm_penalty = c(0.01, 0.1, 1)
  )
)
```

Multivariate series are **not** full tensor products by default. The additive
construction concatenates univariate terms. `combine = "anova"` adds tensor
products only up to `interaction_order`; `combine = "total_degree"` provides
a total-degree Legendre basis. A full `combine = "tensor"` basis is opt-in and
guarded by `max_features` and `max_feature_elements`, because \(J\) terms in
each of \(p\) coordinates can produce order \(J^p\) features. A total-degree
Legendre basis grows more slowly, approximately as \(\binom{p+d}{d}\).

Fold-specific bounds, spline knots, centering constants, and normalization
factors are learned only from the relevant training fold and reused unchanged
for validation and prediction. Fourier bases impose periodic endpoint
behavior. The built-in wavelet is Haar; it is not a
Cohen--Daubechies--Vial boundary-corrected construction. Specialist series can
be supplied as user-defined positive-semidefinite kernels. Such callbacks use
the general kernel route; unlike the built-in series families, they do not
automatically receive the leakage-safe primal feature-map machinery.

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

### Neural quick-start presets

The historical `learner = "relu"` name now denotes the generic neural learner.
Presets provide small candidate grids; explicit detailed controls always take
precedence over preset values.

| `relu_preset` | candidates | activation/backend |
|---|---|---|
| `"affine"` | affine only | native linear |
| `"fast"` | affine, one layer of width 8 | native ReLU |
| `"standard"` | affine, widths 8 and 16 | native ReLU |
| `"flexible"` | affine, widths 8/16, and 16-by-8 | native ReLU and tanh |
| `"nnet"` | affine, widths 4/8/16 | `nnet`, sigmoid, skip connections |

```r
fast_fit <- odrl(
  x, treatment, outcome,
  learner = "relu", loss = "logistic",
  control = odrl_control(relu_preset = "fast")
)

# A preset plus an explicit override:
customized_fit <- odrl(
  x, treatment, outcome,
  learner = "relu", loss = "logistic",
  control = odrl_control(
    relu_preset = "flexible",
    relu_architectures = list(integer(), 12L, c(16L, 8L))
  )
)
```

The native engine accepts an affine candidate (`integer()`), one layer (`8L`),
or arbitrary width vectors such as `c(16L, 8L)`, with ReLU, leaky ReLU, tanh,
sigmoid, or linear activation. These presets are pragmatic tuning grids rather
than claims of universal optimality; policy-value CV (with an optional one-SE
rule) selects among them.

`relu_backend = "nnet"` uses `nnet::nnet()` for an affine model or one sigmoid
hidden layer with optional input-to-output skip connections. It supports the
weighted **logistic** margin loss only; it does not faithfully optimize hinge,
exponential, squared-hinge, or arbitrary custom losses. An affine `nnet`
candidate requires `skip = TRUE`, and wide fits may require a larger `MaxNWts`.
The BFGS-based `nnet` engine is not a deep/GPU backend, and its decay scale is
backend-specific. Use the native engine or named custom `fit`/`predict`
callbacks for other objectives or architectures. Likewise, `tree_backend` may
provide callbacks for another tree implementation.

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
