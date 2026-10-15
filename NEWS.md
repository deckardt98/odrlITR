# odrlITR (development version)

* Made the fitting and tuning code easier to follow with descriptive names,
  explicit statistical steps, and shared helpers for repeated logic. Kept
  public interfaces, defaults, and numerical results unchanged.

* Rewrote the README and getting started guide around a complete workflow,
  with examples for matrices, data frames, named treatments, and prediction.
* Explained treatment recommendations, decision scores, double residual
  scores, nuisance cross-fitting, and policy tuning in plain language.
* Added a source reading guide and development commands for contributors.

# odrlITR 0.1.1

* Updated the hinge SVM with bounded scores to fit the penalized hinge problem
  before hard tanh clipping, matching Example 3 for Gaussian, linear, and
  other supported kernels.
* Simplified the README, vignette, and function documentation and added a
  runnable installation example for the public GitHub package.
* Added the arXiv paper to the package metadata, citation, README, vignette,
  and package help.

# odrlITR 0.1.0

* Initial package release.
* Added configurable, treatment-stratified cross-fitting for Super Learner and
  included parametric nuisance models, plus fold assignments and nuisance
  predictions supplied by users.
* Added direct optimization of linear rules with auditable HiGHS stopping
  controls and optimization of trees with configurable depth through
  `policytree` or custom backends.
* Added hinge, exponential, logistic, squared hinge, and custom differentiable
  surrogate losses for linear, Gaussian, polynomial, and custom RKHS kernels
  and neural network policies.
* Added Legendre, Fourier, B-spline, Haar wavelet, and partitioned local
  polynomial series constructed separately within each fold, with additive,
  low-order ANOVA, total-degree, and optional full tensor-product constructions.
  These use a primal fit rather than a dense kernel matrix for training.
* Added affine, shallow, and multilayer neural policies with common activation
  functions, small architecture presets, tuning by the policy criterion,
  custom losses, an optional `nnet` backend with one hidden layer, and external
  backend callbacks.
