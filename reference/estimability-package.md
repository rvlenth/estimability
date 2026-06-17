# Estimability Tools for Linear Models

Provides tools for determining estimability of linear functions of
regression coefficients, and alternative `epredict` methods for `lm`,
`glm`, and `mlm` objects that handle non-estimable cases correctly.

## Details

|          |                      |
|----------|----------------------|
| Package: | estimability         |
| Type:    | Package              |
| Details: | See DESCRIPTION file |

When a linear model is not of full rank, the regression coefficients are
not uniquely estimable. However, the predicted values are unique, as are
other linear combinations where the coefficients lie in the row space of
the data matrix. Thus, estimability of a linear function of regression
coefficients can be determined by testing whether the coefficients lie
in this row space – or equivalently, are orthogonal to the corresponding
null space.

This package provides functions [`nonest.basis`](nonest.basis.md) and
[`is.estble`](nonest.basis.md) to facilitate such an estimability test.
Package developers may find these useful for incorporating in their
`predict` methods when new predictor settings are involved.

The function [`estble.subspace`](estble.subspace.md) is useful for
projecting a matrix onto an estimable subspace whose rows are all
estimable.

The package also provides [`epredict`](epredict.md) methods –
alternatives to the [`predict`](https://rdrr.io/r/stats/predict.html)
methods in the stats package for `"lm"`, `"glm"`, and `"mlm"` objects.
When the `newdata` argument is specified, estimability of each new
prediction is checked and any non-estimable cases are replaced by `NA`.

## Author

Russell V. Lenth \<russell-lenth@uiowa.edu\>

## References

Monahan, John F. (2008) *A Primer on Linear Models*, CRC Press. (Chapter
3)
