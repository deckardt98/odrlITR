# odrlITR 0.1.0

* Initial package release.
* Added configurable, treatment-stratified cross-fitting for Super Learner and
  built-in parametric nuisance models, plus user-specified fold assignments and
  user-supplied nuisance predictions.
* Added direct linear-rule optimization with auditable HiGHS stopping controls
  and shallow-tree optimization through `policytree` or custom backends.
* Added hinge, logistic, squared-hinge, and custom differentiable surrogate
  losses for linear, Gaussian, polynomial, and custom RKHS kernels.
* Added affine, shallow, and multilayer neural policies with common activation
  functions, criterion-based tuning, custom losses, and external backend
  callbacks.
