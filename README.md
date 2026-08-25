# odrlITR

<!-- badges: start -->
[![R-CMD-check](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`odrlITR` learns individualized treatment rules with orthogonal double
residual learning. It:

1. estimates the propensity score and marginal outcome regression with
   cross-fitting, or accepts user-supplied nuisance predictions;
2. fits a treatment rule from the score
   `(A - e_hat(X)) * (Y - m_hat(X))`.

Treatment is represented internally as `-1` and `+1`; larger outcomes are
preferred.

## Installation

Install the public GitHub package with:

```r
install.packages("pak")
pak::pak("deckardt98/odrlITR")
```

The core SVM and parametric nuisance learners use base R. Other learners use
optional packages:

```r
install.packages(c("SuperLearner", "glmnet", "policytree", "highs", "nnet"))
```

## Quick start

This example has no optional package dependencies:

```r
library(odrlITR)

dat <- odrl_simulate(n = 250, boundary = "linear", seed = 1)

fit <- odrl(
  dat$x, dat$a, dat$y,
  learner = "svm",
  loss = "hinge",
  nuisance = "parametric",
  nuisance_folds = 5,
  control = odrl_control(seed = 1)
)

head(predict(fit, dat$x))
range(predict(fit, dat$x, type = "score"))
summary(fit)
```

For factor or character treatment, set `positive` to the level that should be
treated as `+1`.

## Nuisance estimation

The default first stage is cross-fitted Super Learner:

```r
fit <- odrl(
  dat$x, dat$a, dat$y,
  learner = "svm",
  nuisance_folds = 5,
  sl.library = c("SL.mean", "SL.glm", "SL.glmnet")
)
```

Use `nuisance = "parametric"` for cross-fitted logistic propensity and linear
outcome models. In a randomized trial, pass `known_pi` or `known_e` to use the
known treatment mechanism while cross-fitting the outcome model.

Externally estimated predictions can be supplied with
`odrl_nuisance_user()`:

```r
nuisance <- odrl_nuisance_user(
  m = m_oof,
  pi = pi_oof,
  fold_id = fold_id,
  source = "external cross-fitting",
  out_of_fold = TRUE
)

fit <- odrl(dat$x, dat$a, dat$y, learner = "svm", nuisance = nuisance)
```

## Policy or ITR learners

| `learner` | allowed `loss` | policy class | optional package |
|---|---|---|---|
| `"tree"` | `"exact"` | shallow decision tree using `policytree` or a custom backend | `policytree` |
| `"linear"` | `"exact"` | linear decision function via mixed-integer programming | `highs` |
| `"svm"` | hinge, exponential, logistic, squared hinge, or custom | linear, Gaussian, polynomial, wavelet, Fourier, or custom RKHS score class | none |
| `"relu"` | hinge, exponential, logistic, squared hinge, or custom | configurable feed-forward network or external neural backend | `nnet` for its optional logistic-only backend |

The tree and linear learners directly optimize the empirical double residual
objective, which is equivalent to a cost-sensitive classification risk. SVM
and neural learners fit surrogate losses; their tuning parameters are selected
by cross-validation using the same double residual objective on held-out folds.

## Bounded hinge SVM

For the default Gaussian kernel with hinge loss, `odrlITR` minimizes the
penalized hinge-loss surrogate of the empirical double residual objective,
then applies hard-tanh clipping so that the fitted score lies in `[-1, 1]`:

```text
T1(f) = max(-1, min(1, f))
```

This is the construction in Example 3 of the paper. Clipping preserves the
treatment rule and cannot increase hinge loss. The clipped score is bounded
but need not remain an RKHS function. The paper uses this bounded score for
its universal Neyman-orthogonality result.

Set `svm_hinge_mode = "bounded"` to apply the same fit-then-clip construction
to another supported kernel, including the linear kernel. Set
`svm_hinge_mode = "regularized"` to return the ordinary, unclipped hinge
score.

```r
linear_fit <- odrl(
  dat$x, dat$a, dat$y,
  learner = "svm", loss = "hinge",
  nuisance = "parametric",
  control = odrl_control(
    svm_kernel = "linear",
    svm_hinge_mode = "bounded",
    svm_penalty = c(0.01, 0.1, 1),
    seed = 2
  )
)
```

The fitted score includes a penalized constant term, equivalently using the
augmented kernel `K + 1`. The selected penalty, convergence result, norm, and
primal-dual gap are available in `fit$policy$selected` and
`fit$policy$diagnostics`.

## Other kernels and finite series

Set `svm_kernel` to `"linear"`, `"rbf"`, or `"polynomial"`, or supply a
custom positive-semidefinite kernel. Built-in finite-series score classes are
available for Legendre, Fourier, B-spline, Haar-wavelet, and local-polynomial
bases:

```r
series <- odrl_series_kernel(
  basis = "legendre",
  legendre_degree = c(1, 2, 3),
  combine = "anova",
  interaction_order = 2
)

series_fit <- odrl(
  dat$x, dat$a, dat$y,
  learner = "svm", loss = "logistic",
  nuisance = "parametric",
  control = odrl_control(svm_kernel = series)
)
```

Gaussian and custom-kernel fits retain a dense training kernel matrix.
Finite-series learners use an explicit feature map and avoid that matrix.

## Neural policies

`learner = "relu"` supports affine, shallow, and multilayer scores with ReLU,
leaky-ReLU, tanh, sigmoid, or linear activations. Presets provide compact
candidate grids:

```r
neural_fit <- odrl(
  dat$x, dat$a, dat$y,
  learner = "relu", loss = "logistic",
  nuisance = "parametric",
  control = odrl_control(relu_preset = "fast", seed = 3)
)
```

Available presets are `"affine"`, `"fast"`, `"standard"`, `"flexible"`,
and `"nnet"`. Explicit architecture and optimizer settings override a preset.

## Fitted objects

Use `predict(fit)` for treatment recommendations and
`predict(fit, type = "score")` for decision scores. `summary(fit)` reports the
nuisance and policy fit; learner-specific tuning and optimization details are
stored in `fit$policy$diagnostics`.

See `vignette("getting-started", package = "odrlITR")` and the
[function reference](https://github.com/deckardt98/odrlITR/tree/main/man) for
the complete interface.

`odrlITR` is research software under active development.
