# Find an estimable subspace

Determine a transformation `B` of the rows of a matrix `L` such that
`B %*% L` is estimable. A practical example is in jointly testing a set
of contrasts `L` in a linear model, and we need to restrict to the
subspace spanned by the rows of `L` that are estimable.

## Usage

``` r
estble.subspace(L, nbasis, tol = 1e-08)
```

## Arguments

- L:

  A matrix of dimensions *k* by *p*

- nbasis:

  A *k* by *b* matrix whose columns form a basis for non-estimable
  linear functions – such as is returned by
  [`nonest.basis`](nonest.basis.md)

- tol:

  Numeric tolerance for assessing nonestimability. See
  [`is.estble`](nonest.basis.md).

## Value

An *r* by *p* matrix `M = B %*% L` whose rows are all orthogonal to the
columns of `nbasis`. The matrix `B` is attached as `attr(M, "B")`. Note
that if any rows of `L` were non-estimable, then *r* will be less than
*k*. In fact, if there are no estimable functions in the row space of
`L`, then *r* = 0.

## Details

We require `B` such that all the rows of `M = B %*% L` are estimable,
i.e. orthogonal to the columns of `nbasis`. Thus, we need
`B %*% L %*% nbasis` to be zero, or equivalently, `t(B)` must be in the
null space of `t(L %*% nbasis)`. This can be found using
[`nonest.basis`](nonest.basis.md).

## Examples

``` r
### Find a set of estimable interaction contrasts for a 3 x 4 design 
### with two empty cells.
des <- expand.grid(A = factor(1:3), B = factor(1:4))
des <- des[-c(5, 12), ]  # cells (2,2) and (3,4) are empty

X <- model.matrix(~ A * B, data = des)
N <- nonest.basis(X)

L <- cbind(matrix(0, nrow = 6, ncol = 6), diag(6))
# i.e., give nonzero weight only to interaction effects

estble.subspace(L, N)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]    0    0    0    0    0    0    0    0    1     0     0     0
#> [2,]    0    0    0    0    0    0    0    0    0     1     0     0
#> [3,]    0    0    0    0    0    0    0    0    0     0     1     0
#> [4,]    0    0    0    0    0    0    0   -1    0     0     0     0
#> attr(,"B")
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    0    1    0    0    0
#> [2,]    0    0    0    1    0    0
#> [3,]    0    0    0    0    1    0
#> [4,]    0   -1    0    0    0    0

# Tougher demo: create a variation where all rows of L are non-estimable
set.seed(1921.0923)  # Franklin Graybill's birthday
LL <- matrix(rnorm(36), ncol = 6) %*% L
estble.subspace(LL, N)
#>      [,1] [,2] [,3] [,4] [,5] [,6]          [,7]       [,8]       [,9]
#> [1,]    0    0    0    0    0    0 -7.024565e-17  1.0316251  0.2402223
#> [2,]    0    0    0    0    0    0 -1.914267e-16  0.9679828  0.3871052
#> [3,]    0    0    0    0    0    0  3.613173e-16 -0.8807992 -0.4495615
#> [4,]    0    0    0    0    0    0  4.003655e-16  0.1534469 -1.8498715
#>            [,10]       [,11]         [,12]
#> [1,]  1.12000141 -0.67713307  5.497174e-17
#> [2,] -0.04158445 -1.04336175  9.151698e-17
#> [3,] -1.22043728  0.89366353 -2.091758e-16
#> [4,]  0.63353316  0.01685268 -1.332937e-17
#> attr(,"B")
#>            [,1]        [,2]        [,3]       [,4]       [,5]        [,6]
#> [1,]  0.2368168  0.44441284  0.81010361  0.2315774  0.1846316 -0.04929783
#> [2,] -0.7367046 -0.04759893  0.17253661  0.5347192 -0.3145821  0.20086111
#> [3,] -0.4882847 -0.14729815  0.15061904 -0.3305079  0.7689707  0.12901096
#> [4,]  0.3357518 -0.18734622 -0.01314464  0.1720998  0.1026899  0.90101930
```
