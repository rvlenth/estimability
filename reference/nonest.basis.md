# Estimability Tools

This documents the functions needed to test estimability of linear
functions of regression coefficients.

## Usage

``` r
nonest.basis(x, ...)

# Default S3 method
nonest.basis(x, ...)

# S3 method for class 'qr'
nonest.basis(x, ...)

# S3 method for class 'matrix'
nonest.basis(
  x,
  tol = 5e-08,
  rank = attr(x, "rank"),
  pivot = attr(x, "pivot"),
  ...
)

# S3 method for class 'lm'
nonest.basis(x, ...)

# S3 method for class 'svd'
nonest.basis(x, tol = 5e-08, rank = NULL, pivot = NULL, ...)

legacy.nonest.basis(x, ...)

is.estble(x, nbasis, tol = 1e-08)

all.estble
```

## Arguments

- x:

  For `nonest.basis`, an object of a class in `methods("nonest.basis")`.
  Or, in `is.estble`, a numeric vector or matrix for assessing
  estimability of `sum(x * beta)`, where `beta` is the vector of
  regression coefficients.

- ...:

  Additional arguments passed to other methods.

- tol:

  Numeric tolerance for assessing rank or nonestimability. For
  determining rank, singular values less than `tol` times the largest
  singular value are regarded as zero. For determining estimability with
  a nonzero \\x\\, \\\beta'x\\ is assessed by whether or not
  \\\|\|N'x\|\|^2 \< \tau \|\|x'x\|\|^2\\, where \\N\\ and \\\tau\\
  denote `nbasis` and `tol`, respectively.

- rank:

  Optional integer value. If not `NULL`, it is used in place of `tol` to
  decide the rank of the model matrix

- pivot:

  Optional integer vector of length equal to the number of columns in
  the model matrix. If non-`NULL`, it is assumed that this specifies the
  order of columns in `x` relative to the order of predictors in the
  model matrix, and the returned non-estimable basis is permuted
  accordingly.

- nbasis:

  Matrix whose columns span the null space of the model matrix. Such a
  matrix is returned by `nonest.basis`.

## Value

When \\X\\ is not full-rank, the methods for `nonest.basis` return a
basis for the null space of \\X\\. The number of rows is equal to the
number of regression coefficients (*including* any `NA`s); and the
number of columns is equal to the rank deficiency of the model matrix.
The columns are orthonormal. If the model is full-rank, then
`nonest.basis` returns `all.estble`. The `matrix` method uses \\X\\
itself, the `qr` method uses the \\QR\\ decomposition of \\X\\, and the
`lm` method recovers the required information from the object.

The function `is.estble` returns a logical value (or vector, if `x` is a
matrix) that is `TRUE` if the function is estimable and `FALSE` if not.

## Details

Consider a linear model \\y = X\beta + E\\. If \\X\\ is not of full
rank, it is not possible to estimate \\\beta\\ uniquely. However,
\\X\beta\\ *is* uniquely estimable, and so is \\a'X\beta\\ for any
conformable vector \\a\\. Since \\a'X\\ comprises a linear combination
of the rows of \\X\\, it follows that we can estimate any linear
function where the coefficients lie in the row space of \\X\\.
Equivalently, we can check to ensure that the coefficients are
orthogonal to the null space of \\X\\.

The `nonest.basis` method for class `'svd'` is not really functional as
a method because there is no `"svd"` class (at least in R \<= 4.2.0).
But the function `nonest.basis.svd` is exported and may be called
directly; it works with results of
[`svd`](https://rdrr.io/r/base/svd.html) or
[`La.svd`](https://rdrr.io/r/base/svd.html). We *require* `x$v` to be
the complete matrix of right singular values; but we do not need `x$u`
at all.

The `default` method does serve as an `svd` method, in that it only
works if `x` has the required elements of an SVD result, in which case
it passes it to `nonest.basis.svd`.

The `matrix` method runs `nonest.basis.svd(svd(x, nu = 0))`. The `lm`
method runs the `qr` method on `x$qr`.

The function `legacy.nonest.basis` is the original default method in
early versions of the estimability package. It may be called with `x`
being either a matrix or a `qr` object, and after obtaining the `R`
matrix, it uses an additional QR decomposition of `t(R)` to obtain the
needed basis. (The current `nonest.basis` method for `qr` objects is
instead based on the singular-value decomposition of R, and requires
much simpler code.)

The constant `all.estble` is simply a 1 x 1 matrix of `NA`. This
specifies a trivial non-estimability basis, and using it as `nbasis`
will cause everything to test as estimable.

## Use with Cholesky decompositions

If the [`chol`](https://rdrr.io/r/base/chol.html) function is called
with `pivot = TRUE`, the returned matrix is suitable to provide to the
`matrix` method of `nonest.basis`, and the returned basis is correctly
pivoted without specifying `pivot` explicitly.

This works because the Cholesky decomposition \\R\\ of \\X'X\\ (as
constructed by [`chol()`](https://rdrr.io/r/base/chol.html)) satisfies
\\X'X = R'R\\, implying that the rows of \\X'X\\ are linear combinations
of the rows of \\R\\; and they are both of equal rank - so \\R\\ and
\\X'X\\ have the same row space.

## Examples

``` r
require(estimability)

X <- cbind(rep(1,5), 1:5, 5:1, 2:6)
( nb <- nonest.basis(X) )
#>            [,1]        [,2]
#> [1,]  0.4171348 -0.89010713
#> [2,] -0.7077494 -0.26855996
#> [3,] -0.1606977  0.08879245
#> [4,]  0.5470517  0.35735241

# via the singular-value decomposition
SVD <- svd(X, nu = 0)    # we don't need the U part of UDV'
nonest.basis.svd(SVD)    # same result as above
#>            [,1]        [,2]
#> [1,]  0.4171348 -0.89010713
#> [2,] -0.7077494 -0.26855996
#> [3,] -0.1606977  0.08879245
#> [4,]  0.5470517  0.35735241

# Via the Cholesky decomposition of X'X
ch <- chol(t(X) %*% X, pivot = TRUE) |> suppressWarnings()
( nb.chol <- nonest.basis(ch) )
#>            [,1]        [,2]
#> [1,]  0.9733285 -0.13756349
#> [2,] -0.1622214 -0.73940376
#> [3,] -0.1622214 -0.08597718
#> [4,]  0.0000000  0.65342658

# Test estimability of some linear functions for this X matrix
lfs <- rbind(c(1,4,2,5), c(2,3,9,5), c(1,2,2,1), c(0,1,-1,1))
is.estble(lfs, nb)
#> [1]  TRUE  TRUE FALSE  TRUE
is.estble(lfs, nb.chol) # same results even though nb.chol is not identical to nb
#> [1]  TRUE  TRUE FALSE  TRUE

# Illustration on 'lm' objects:
warp.lm1 <- lm(breaks ~ wool * tension, data = warpbreaks,
    subset = -(26:38),
    contrasts = list(wool = "contr.treatment", tension = "contr.treatment"))
zapsmall(warp.nb1 <- nonest.basis(warp.lm1))
#>            [,1]
#> [1,]  0.0000000
#> [2,] -0.5773503
#> [3,]  0.0000000
#> [4,]  0.0000000
#> [5,]  0.5773503
#> [6,]  0.5773503

# different parameterization of the same model:
warp.lm2 <- update(warp.lm1,
    contrasts = list(wool = "contr.sum", tension = "contr.helmert"))
zapsmall(warp.nb2 <- nonest.basis(warp.lm2))
#>            [,1]
#> [1,] -0.3779645
#> [2,]  0.3779645
#> [3,]  0.5669467
#> [4,]  0.1889822
#> [5,] -0.5669467
#> [6,] -0.1889822

# These bases depend of the parameterizations, but they both 
# correctly identify the empty cell
wcells = with(warpbreaks, expand.grid(wool = levels(wool), tension = levels(tension)))
epredict(warp.lm1, newdata = wcells, nbasis = warp.nb1)
#>        1        2        3        4        5        6 
#> 44.55556       NA 24.00000 27.28571 25.71429 18.77778 
epredict(warp.lm2, newdata = wcells, nbasis = warp.nb2)
#>        1        2        3        4        5        6 
#> 44.55556       NA 24.00000 27.28571 25.71429 18.77778 
```
