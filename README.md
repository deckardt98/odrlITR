# odrl

<!-- badges: start -->
[![R-CMD-check](https://github.com/deckardt98/odrl/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/deckardt98/odrl/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`odrl` implements **orthogonal double residual learning** for binary
individualized treatment rules. The name is short, searchable, and matches the
method used in the paper without tying the software to one optimizer or policy
class.

The package separates two stages:

1. estimate `pi(x) = P(A = +1 | X = x)` and `m(x) = E(Y | X = x)` by
   five-fold cross-fitted Super Learner, or supply aligned nuisance predictions;
2. optimize a policy using the signed score
   `(A - e_hat(X)) * (Y - m_hat(X))`, where `e_hat = 2 * pi_hat - 1`.

## Installation

```r
# Once the repository is public:
# install.packages("pak")
# pak::pak("deckardt98/odrl")

# From a downloaded CRAN-style source archive:
install.packages("odrl_0.1.0.tar.gz", repos = NULL, type = "source")
```

Optional engines are installed only when needed:

```r
install.packages(c("SuperLearner", "glmnet", "policytree", "highs"))
```

## Quick start

### Built-in cross-fitted Super Learner

```r
library(odrl)

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

Propensity and outcome libraries can be chosen separately with
`sl.library.pi` and `sl.library.m`. In a randomized trial, pass `known_pi`
(or `known_e`) to retain the known assignment mechanism while still
cross-fitting `m`.

For factor or character treatment, use `positive = "treated"` to make the
`+1` action explicit.

### User-supplied nuisances

```r
nuisance <- odrl_nuisance_user(
  m = m_oof,
  pi = pi_oof,
  fold_id = fold_id,
  source = "external five-fold cross-fitting",
  out_of_fold = TRUE
)

fit <- odrl(x, treatment, outcome, learner = "svm", loss = "logistic",
            nuisance = nuisance)
```

The nuisance argument may also be a list or a function returning a nuisance
object. No propensity clipping is applied by default because ODRL does not
divide by the propensity.

## Policy learners

| `learner` | allowed `loss` | policy class | optional package |
|---|---|---|---|
| `"tree"` | `"exact"` | shallow axis-aligned tree over `policytree`'s candidate splits | `policytree` |
| `"linear"` | `"exact"` | bounded-margin affine sign rule via mixed-integer optimization | `highs` |
| `"svm"` | `"hinge"`, `"logistic"` | Gaussian RKHS (or linear kernel for logistic) | none |
| `"relu"` | `"hinge"`, `"logistic"` | affine or one-hidden-layer ReLU score | none |

The bounded-hinge SVM uses a normalized Gaussian kernel and an RKHS radius at
most one, which certifies scores in `[-1, 1]`. The ReLU hinge learner applies a
hard-tanh output map. Logistic loss is classification-calibrated, but it does
not have the bounded-hinge objective's special universal orthogonality. Its
constant offset and RKHS component are both included in the quadratic
penalty.

Gaussian-kernel learners form and retain an \(n \times n\) training kernel
matrix. Their time and memory costs therefore grow at least quadratically in
the training sample size, and prediction requires kernel evaluations against
all training observations. For large samples, prefer the affine, tree, or
ReLU learner, or benchmark the kernel learner on a representative subset
before committing to a full tuning grid.

Affine optimization is genuinely mixed-integer. A time-limited feasible
incumbent is not described as proven optimal: solver status, MIP gap,
constraint checks, and the objective audit are retained in
`fit$policy$diagnostics`. Set `linear_require_gap = TRUE` when a certificate is
required.

`linear_require_gap = TRUE` certifies only the requested relative gap. For a
numerical global-optimum certificate within the bounded-margin standardized
affine class, also set `linear_relative_gap = 0`; inspect
`certified_global_optimum` in the diagnostics.

For numerical conditioning, policy optimization and tuning use the score
divided by its sample mean absolute value. This positive rescaling leaves the
unpenalized binary objective and treatment rule unchanged and makes the
package's penalty grids interpretable on a common scale. The reported training
criterion is recomputed from the original, unscaled score.

```r
control <- odrl_control(
  tree_depth = 2,
  linear_time_limit = 300,
  linear_relative_gap = 0.01,
  svm_rbf_multiplier = c(0.5, 1, 2),
  svm_penalty = c(0.01, 0.1, 1),
  relu_hidden_units = c(0, 8, 16),
  relu_decay = c(0.001, 0.01, 0.1),
  relu_selection = "one_se",
  seed = 2026
)
```

## Reproducibility and scope

- Outer nuisance cross-fitting and second-stage policy tuning use separate
  folds.
- Treatment ties and zero scores use the paper's `+1` convention.
- Raw and numerically thresholded scores are both retained for auditing.
- The package learns a policy; causal interpretation still requires the usual
  identification conditions and appropriate nuisance estimation.
- The ReLU engine is a regularized one-hidden-layer subclass with an optional
  affine candidate. Its hard-tanh hinge output has the paper's bounded-score
  orthogonality property, but the package does not claim the manuscript's
  sparse deep-network rate for this particular optimizer.

See `vignette("getting-started", package = "odrl")` for a complete workflow.

## Development status

Version 0.1.0 is an initial research-software release prepared for public
review and eventual CRAN submission. Please report bugs through the GitHub
issue tracker.
