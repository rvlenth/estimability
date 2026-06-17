# Estimability Enhancements for `lm` and Relatives

These functions call the corresponding S3 `predict` methods in the stats
package, but with a check for estimability of new predictions, and with
appropriate actions for non-estimable cases.

## Usage

``` r
epredict(object, ...)

# S3 method for class 'lm'
epredict(
  object,
  newdata,
  ...,
  type = c("response", "terms", "matrix", "estimability"),
  nonest.tol = 1e-08,
  nbasis = object$nonest
)

# S3 method for class 'glm'
epredict(
  object,
  newdata,
  ...,
  type = c("link", "response", "terms", "matrix", "estimability"),
  nonest.tol = 1e-08,
  nbasis = object$nonest
)

# S3 method for class 'mlm'
epredict(
  object,
  newdata,
  ...,
  type = c("response", "matrix", "estimability"),
  nonest.tol = 1e-08,
  nbasis = object$nonest
)

eupdate(object, ...)

# S3 method for class 'lm'
eupdate(object, ...)
```

## Arguments

- object:

  An object inheriting from `lm`

- ...:

  Arguments passed to [`predict`](https://rdrr.io/r/stats/predict.html)
  or [`update`](https://rdrr.io/r/stats/update.html)

- newdata:

  A `data.frame` containing predictor combinations for new predictions

- type:

  haracter string specifying the desired result. See Details.

- nonest.tol:

  Tolerance used by [`is.estble`](nonest.basis.md) to check estimability
  of new predictions

- nbasis:

  a basis for the null space, e.g., a result of a call to
  [`nonest.basis`](nonest.basis.md). If `nbasis` is `NULL`, a basis is
  constructed from `object`.

## Value

The same as the result of a call to the `predict` method in the stats
package, except rows or elements corresponding to non-estimable
predictor combinations are set to `NA`. The value for `type` is
`"matrix"` or `"estimability"` is explained under details.

## Details

If `newdata` is missing or `object` is not rank-deficient, this method
irectly to the same method in the stats library. In rank-deficient cases
with `newdata` provided, each row of `newdata` is tested for
estimability against the null basis provided in `nbasis`. Any
non-estimable cases found are replaced with `NA`s.

The `type` argument is passed to
[`predict`](https://rdrr.io/r/stats/predict.html) when it is one of
`"response"`, `"link"`, or `"terms"`. With `newdata` present and
`type = "matrix"`, the model matrix for `newdata` is returned, with an
attribute `"estble"` that is a logical vector of length `nrow(newdata)`
indicating whether each row is estimable. With `type = "estimability"`,
just the logical vector is returned.

If you anticipate making several `epredict` calls with new data, it
improves efficiency to either obtain the null basis and provide it in
the call, or add it to `object` with the name `"nonest"` (perhaps via a
call to `eupdate`).

`eupdate` is an S3 generic function with a method provided for `"lm"` It
updates the object according to any arguments in `...`, then obtains the
updated object's nonestimable basis and returns it in `object$nonest`.

## Note

The capabilities of the `epredict` function for `lm` objects is provided
by [`predict.lm`](https://rdrr.io/r/stats/predict.lm.html) (if using R
version 4.3.0 or later) with `rankdeficient = "NA"`; however, `epredict`
uses estimability's own criteria to determine which predictions are set
to `NA`. An advantage of using `epredict` is one of efficiency: we can
compute the null basis once and for all and have it available additional
predictions, whereas `predict.lm` will re-compute it each time. If the
user wishes to see a message explaining why `NA`s were displayed, set
`options(estimability.verbose = TRUE)`.

## Examples

``` r
require("estimability")

# Fake data where x3 and x4 depend on x1, x2, and intercept
x1 <- -4:4
x2 <- c(-2,1,-1,2,0,2,-1,1,-2)
x3 <- 3*x1 - 2*x2
x4 <- x2 - x1 + 4
y <- 1 + x1 + x2 + x3 + x4 + c(-.5,.5,.5,-.5,0,.5,-.5,-.5,.5)

# Different orderings of predictors produce different solutions
mod1234 <- lm(y ~ x1 + x2 + x3 + x4)
mod4321 <- eupdate(lm(y ~ x4 + x3 + x2 + x1))
# (Estimability checking with mod4321 will be more efficient because
#  it will not need to recreate the basis)
mod4321$nonest
#>             [,1]        [,2]
#> [1,] -0.17177076  0.94865978
#> [2,]  0.04294269 -0.23716495
#> [3,] -0.24764833 -0.13231961
#> [4,] -0.53823935 -0.02747428
#> [5,]  0.78588768  0.15979389

 
# test data:
testset <- data.frame( 
              x1 = c(3,  6,  6,  0,  0,  1), 
              x2 = c(1,  2,  2,  0,  0,  2), 
              x3 = c(7, 14, 14,  0,  0,  3), 
              x4 = c(2,  4,  0,  4,  0,  4))

# Look at predictions when we don't check estimability
suppressWarnings( # Disable the warning from stats::predict.lm
    rbind(p1234 = predict(mod1234, newdata = testset),
          p4321 = predict(mod4321, newdata = testset)))
#>        1  2  3 4   5  6
#> p1234 14 23 23 5   5  8
#> p4321 14 47 23 5 -19 14

# Compare with results when we do check:
rbind(p1234 = epredict(mod1234, newdata = testset),
      p4321 = epredict(mod4321, newdata = testset))
#>        1  2  3 4  5  6
#> p1234 14 NA 23 5 NA NA
#> p4321 14 NA 23 5 NA NA

# stats::predict has same capability for lm objects starting in version 4.3.0:
if((R.Version()$major >= 4) && (R.Version()$minor >= 3))
  stats::predict(mod1234, newdata = testset, rankdeficient = "NA")
#>  1  2  3  4  5  6 
#> 14 NA 23  5 NA NA 

# Note that estimable cases have the same predictions

# change mod1234 and include nonest basis 
mod134 <- eupdate(mod1234, . ~ . - x2, subset = -c(3, 7))
mod134$nonest
#>            [,1]
#> [1,]  0.9561829
#> [2,]  0.1195229
#> [3,] -0.1195229
#> [4,] -0.2390457

# When row spaces are the same, bases are interchangeable
# so long as you account for the ordering of parameters:
epredict(mod4321, newdata = testset, type = "estimability",
    nbasis = nonest.basis(mod1234)[c(1,5:2), ])
#>     1     2     3     4     5     6 
#>  TRUE FALSE  TRUE  TRUE FALSE FALSE 
    
# Comparison with predict.lm in R >= 4.3.0
if((R.Version()$major >= 4) && (R.Version()$minor >= 3))
  stats::predict(mod4321, newdata = testset, rankdeficient = "NA")
#>  1  2  3  4  5  6 
#> 14 NA 23  5 NA NA 

if (FALSE) { # \dontrun{
### Additional illustration
example(nonest.basis)  ## creates model objects warp.lm1 and warp.lm2

# The two models have different contrast specs. But the empty cell
# is correctly identified in both:
fac.cmb <- expand.grid(wool = c("A", "B"), tension = c("L", "M", "H"))
cbind(fac.cmb, 
      pred1 = epredict(warp.lm1, newdata = fac.cmb), 
      pred2 = epredict(warp.lm2, newdata = fac.cmb))
} # }
```
