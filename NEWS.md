# odrlITR 0.1.0

* Initial package release.
* Added configurable, treatment-stratified cross-fitting for Super Learner and
  built-in parametric nuisance models, plus user-specified fold assignments and
  user-supplied nuisance predictions.
* Added direct linear-rule optimization with auditable HiGHS stopping controls
  and configurable-depth tree optimization through `policytree` or custom
  backends.
* Added hinge, exponential, logistic, squared-hinge, and custom differentiable
  surrogate losses for linear, Gaussian, polynomial, and custom RKHS kernels
  and neural-network policies.
* Added fold-safe Legendre, Fourier, B-spline, Haar-wavelet, and partitioned
  local-polynomial series with additive, low-order ANOVA, total-degree, and
  opt-in full-tensor multivariate constructions. These use a primal fit rather
  than a dense training-kernel matrix.
* Added affine, shallow, and multilayer neural policies with common activation
  functions, quick-start architecture presets, criterion-based tuning, custom
  losses, an optional single-hidden-layer `nnet` backend, and external backend
  callbacks.
