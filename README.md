# odrlITR

<!-- badges: start -->
[![R-CMD-check](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/deckardt98/odrlITR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`odrlITR` is an R package for learning which of two treatments to recommend from a person's measured characteristics. You supply the characteristics, the treatment received, and the observed outcome. The package fits an individualized treatment rule that can recommend a treatment for new observations. Larger outcomes are treated as better.

The method is **orthogonal double residual learning (ODRL)**. It first estimates how treatment and outcome vary with the characteristics, then uses the two residuals to learn a treatment rule. You can choose a decision tree, a linear rule, a kernel method, or a neural network.

## Install

```r
install.packages("pak")
pak::pak("deckardt98/odrlITR")
```

The example below uses the included parametric models and SVM learner. It needs no optional modeling packages.

## Fit your first rule

```r
library(odrlITR)

# Generate an example with five characteristics and two treatments.
dat <- odrl_simulate(n = 200, boundary = "linear", seed = 1)

fit <- odrl(
  x = dat$x,
  a = dat$a,
  y = dat$y,
  learner = "svm",
  nuisance = "parametric",
  control = odrl_control(
    svm_penalty = c(0.05, 0.2),
    svm_rbf_multiplier = 1,
    seed = 1
  )
)

# Recommend treatments for new observations.
new_dat <- odrl_simulate(n = 50, boundary = "linear", seed = 2)
head(predict(fit, newdata = new_dat$x))

# Inspect the fit and its selected tuning settings.
summary(fit)
```

The three inputs have one row or entry per observation, in the same order.

| Input | Meaning | Accepted form |
|---|---|---|
| `x` | Characteristics available when choosing treatment | Numeric matrix or data frame, including factor columns |
| `a` | Treatment actually received | `0`/`1`, `-1`/`+1`, logical, or two treatment labels |
| `y` | Observed outcome, with larger values preferred | Numeric vector |

Handle missing values before fitting. For factor or character treatments, use `positive` to name the treatment represented internally by `+1`. Recommendations use the original treatment labels or numeric coding.

`predict(fit, newdata)` returns a treatment recommendation for each new row. `fitted(fit)` returns recommendations for the training rows. For an SVM, `predict(fit, newdata, type = "score")` returns the decision score whose sign determines the recommendation. This score is on the learner's scale, not the outcome or treatment effect scale.

## Use your own data frame

The formula interface names the outcome, treatment column, and predictors. This runnable example converts the simulated data into the same layout you would use for your own data.

```r
trial <- data.frame(
  outcome = dat$y,
  treatment = ifelse(dat$a == 1, "active", "control"),
  dat$x
)

formula_fit <- odrl_formula(
  outcome ~ X1 + X2 + X3 + X4 + X5,
  treatment = "treatment",
  data = trial,
  positive = "active",
  learner = "svm",
  nuisance = "parametric",
  control = odrl_control(
    svm_penalty = 0.1,
    svm_rbf_multiplier = 1,
    seed = 1
  )
)

# Prediction needs only the predictor columns.
head(predict(formula_fit, newdata = as.data.frame(new_dat$x)))
```

Replace `trial` and the column names with your data. Choose predictors that are available before the treatment decision. You can use transformations and interactions in the formula. With `outcome ~ .`, the treatment column is automatically excluded from the predictors.

## Choose the next step

**Choose the rule you want to use.** The `learner` argument controls the form of the treatment rule.

| `learner` | Form of the rule | Extra package |
|---|---|---|
| `"tree"` | A sequence of splits, with a treatment in each final node | `policytree` |
| `"linear"` | The sign of a weighted sum of predictors and an intercept | `highs` |
| `"svm"` | A score built from a kernel or finite series | None |
| `"relu"` | An affine score or neural network | None for the native backend; `nnet` for that backend |

Install an extra package only when using that learner, for example `install.packages("policytree")`. Tree and linear learners optimize the double residual criterion directly and use `loss = "exact"`. SVM and neural learners fit a surrogate loss and select tuning values by cross-validation. Their default loss is `"hinge"`; `"logistic"`, `"exponential"`, and `"squared_hinge"` are also available.

**Choose the first-stage models.** These models are called *nuisance models* because they support learning the treatment rule. `nuisance = "parametric"` uses logistic regression for treatment and linear regression for outcome. The default, `nuisance = NULL`, uses Super Learner and requires `SuperLearner` plus the packages used by its library, including `glmnet` for the default library. Both options use cross-fitting to predict each row from models trained on other rows. You can also supply your own predictions with `odrl_nuisance_user()`.

For a randomized trial, supply the known allocation probability through `known_pi`, such as `known_pi = 0.5` for equal allocation. This probability refers to the treatment coded as `+1`. The outcome model is still cross-fitted.

The default `odrl()` call uses the tree learner and Super Learner. The examples above set both choices explicitly so they run with the package's included models.

**Adjust complexity when needed.** Use `odrl_control()` to set tree depth, SVM kernels and penalties, or neural network size. The [getting started guide](vignettes/getting-started.Rmd) explains these choices, bounded hinge scores, finite series, custom learners, and how to inspect a fitted object.

## Learn more

- Read the [getting started guide](vignettes/getting-started.Rmd) for a step-by-step workflow and the underlying score.
- Use `?odrl`, `?odrl_control`, and `?odrl_nuisance_user` in R, or browse the [function documentation](man).
- Read [CONTRIBUTING.md](CONTRIBUTING.md) for a map of the source code and instructions for making changes.
- Report problems or suggest improvements in [GitHub issues](https://github.com/deckardt98/odrlITR/issues).

The methodology is described in [“Orthogonal double residual learning for optimal individualized treatment rules”](https://arxiv.org/abs/2608.24085) by Jiaqi Tong and Fan Li (2026), an arXiv preprint. Run `citation("odrlITR")` for the package and paper citations.
