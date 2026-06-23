# Changelog

## estimability 2.0.0

- Changing maintainer to Michael Dumelle.
- Converted `documentation/NAMESPACE` to use `roxygen2`.
- Added (implicit) support for `chol(..., pivot = TRUE)` results. This
  is done by checking for “rank” and “pivot” attributes, and using them
  if present
- Added discussion of Cholesky decompositions to documentation.
- Added a complete example to the vignette.

## estimability 1.5.1

CRAN release: 2024-05-12

- Added pkgdown site.
- We rolled back the required version of R to 4.1.0, same as `emmeans`
  package, so as to not require `>=4.3.0` in packages that depend on
  this one.

## estimability 1.5

CRAN release: 2024-02-20

- We now require R \>= 4.3.0, which plays along with the estimability
  changes in [`predict.lm()`](https://rdrr.io/r/stats/predict.lm.html)
  that came with R 4.3.0.
- We re-coded the `nonest.basis.qr` to something much simpler, and the
  old version is kept available as `legacy.nonest.basis`.
- Added a vignette was added to help developers add estimability
  checking to their package.

## estimability 1.4.1

CRAN release: 2022-08-05

- Correction to version 1.4. The new svd-based methods worked correctly
  only for n x p matrices with n \>= p. Otherwise things go badly awry.
  And this is a big problem because I replaced the default
  [`nonest.basis()`](../reference/nonest.basis.md) method with the svd
  method.

## estimability 1.4

CRAN release: 2022-07-03

- Added support for results of
  [`svd()`](https://rdrr.io/r/base/svd.html), via `nonest.basis.svd`
  function and `default` method. The `matrix` method now uses the SVD
  instead of the QR decomposition.

## estimability 1.3

CRAN release: 2018-02-11

- Added [`estble.subspace()`](../reference/estble.subspace.md) function.

## estimability 1.2.1

- Moved codebase to github repository rvlenth/estimability.

## estimability 1.2

CRAN release: 2016-11-19

- Modified license to make it more compatible with dependents.

## estimability 1.1-1

CRAN release: 2015-09-03

- Added imports of non-base packages that are referenced.

## estimability 1.1

CRAN release: 2015-02-12

- Design improvements to aid in potential scope and usability:
  - Made [`nonest.basis()`](../reference/nonest.basis.md) a generic,
    with provided methods for “qr”, “matrix”, and “lm”.
  - Added [`eupdate()`](../reference/epredict.md) generic and “lm”
    method for updating a model object and including its nonestimability
    basis as part of the object. Added `type = "matrix"` and
    `type = "estimability"` options for
    [`epredict()`](../reference/epredict.md).
