%    Copyright (c) 2015-2016 Russell V. Lenth                                #

\name{estimability-package}
\alias{estimability-package}
\alias{estimability}
\docType{package}
\title{
Estimability Tools for Linear Models
}
\description{
Provides tools for determining estimability of linear functions of regression coefficients, 
and alternative \code{epredict} methods for \code{lm}, \code{glm}, and \code{mlm} objects that handle non-estimable cases correctly.
}
\details{
\tabular{ll}{
Package: \tab estimability\cr
Type: \tab Package\cr
Details: \tab See DESCRIPTION file\cr
}
When a linear model is not of full rank, the regression coefficients are not uniquely estimable. However, the predicted values are unique, as are other linear combinations where the coefficients lie in the row space of the data matrix. Thus, estimability of a linear function of regression coefficients can be determined by testing whether the coefficients lie in this row space -- or equivalently, are orthogonal to the corresponding null space. 
This package provides functions \code{\link{nonest.basis}} and \code{\link{is.estble}} to facilitate such an estimability test. Package developers may find these useful for incorporating in their \code{predict} methods when new predictor settings are involved.

The package also provides \code{\link{epredict}} methods -- alternatives to the \code{\link{predict}} methods in the \pkg{stats} package for \code{"lm"}, \code{"glm"}, and \code{"mlm"} objects. When the \code{newdata} argument is specified, estimability of each new prediction is checked and any non-estimable cases are replaced by \code{NA}. 
}
\author{
Russell V. Lenth <russell-lenth@uiowa.edu>
}
\references{
Monahan, John F. (2008) \emph{A Primer on Linear Models}, CRC Press. (Chapter 3)
}

\keyword{ package }
\keyword{ models }
\keyword{ regression }
