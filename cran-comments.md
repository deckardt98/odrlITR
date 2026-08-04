## R CMD check results

0 errors | 0 warnings | 3 notes

## Test environments

* local macOS, R 4.5.2

## Notes

This is a new submission under private repository review. The repository and
pkgdown URLs are intentionally inaccessible to unauthenticated users and will
be made public before a CRAN submission.

The local check process did not discover Pandoc while checking the top-level
Markdown files, although the package vignettes built and rebuilt successfully
using the bundled Quarto Pandoc executable.

HTML validation and math-rendering checks were skipped because a recent HTML
Tidy executable and the optional `V8` package are not installed locally.

Optional optimization and nuisance-learning engines are listed in `Suggests`
and are checked at runtime with informative errors.
