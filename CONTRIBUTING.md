# Contributing to odrlITR

The goal is code that researchers can read, understand, and adapt. Prefer a clear sequence of steps, descriptive names, and small examples. Open an issue to discuss substantial changes so the method and user interface can be considered together.

## Find your way through the code

Start with the example in [README.md](README.md), then follow the data through these files.

| File | What to look for |
|---|---|
| `R/fit.R` | `odrl()` and `odrl_formula()`, from input preparation through the fitted result |
| `R/nuisance.R` | Cross-fitting the treatment and outcome models, or accepting supplied predictions |
| `R/score.R` | Constructing the double residual score and its classification labels and weights |
| `R/exact-tree.R` and `R/exact-linear.R` | Fitting a tree or affine rule directly to the score |
| `R/svm.R` | Building candidate kernel models, comparing them, and refitting the selected model |
| `R/series-kernel.R` | Creating finite series features for the SVM learner |
| `R/relu.R` and `R/nnet-backend.R` | Neural network fitting, tuning, and the optional `nnet` backend |
| `R/methods.R` | Prediction, fitted actions, coefficients, and printed summaries |
| `R/control.R` | Public settings and neural network presets |
| `R/utils.R` | Shared input encoding, standardization, losses, and numerical helpers |
| `R/simulate.R` | A small data generator with known nuisance functions and treatment rule |

The central quantity is `(a_pm1 - e) * (y - m)`. Internally, treatments are `-1` and `+1`. Policy fitting uses scores divided by their mean absolute value; the reported training criterion uses the raw scores. Preserve the treatment mapping, row alignment, and score scale when changing the fitting path.

## Set up a local copy

Clone the repository and open R in the package directory. Install the tools and dependencies needed for development.

```r
install.packages(c("pak", "devtools", "roxygen2"))
pak::local_install_dev_deps()
devtools::load_all()
```

Run an example from the README to see the current behavior before editing. `devtools::load_all()` reloads your local source without installing the package.

## Make a readable change

- Name variables for their role, such as training rows, validation scores, or selected candidate. Keep familiar mathematical notation where it helps connect the code to the method.
- Keep the main fitting flow easy to follow. Extract a helper when it gives a meaningful step a name or removes repeated logic.
- Use comments to explain a mathematical choice, a transformation, or a numerical constraint. Let the code show routine assignments and loops.
- Put user-facing explanations and examples in the function's roxygen comments. Update the README or vignette when the workflow changes.
- Add a focused test when a change affects behavior. Existing tests in `tests/testthat/` are grouped by learner or feature.

For optimization changes, preserve the signed-score objective and inspect convergence and solver diagnostics. For cross-fitting changes, check that each prediction uses the intended held-out rows. For prediction changes, check both numeric and labeled treatments.

## Check the result

Run these commands in R from the package directory.

```r
# Rebuild help files after editing roxygen comments.
roxygen2::roxygenise()

# Run the relevant tests while iterating.
devtools::test(filter = "nuisance-score")

# Run the full suite and package checks before submitting.
devtools::test()
devtools::check()
```

Replace the test filter with the feature you changed. Rerun the README examples and build the vignette with `devtools::build_vignettes()` when changing documentation. Generated files in `man/` should be committed with the roxygen source that produced them.

Describe the behavior your change affects, why it is easier to use or understand, and which checks you ran. Mention any optional package that prevented a check from running.
