%    Copyright (c) 2015-2016 Russell V. Lenth

\name{nonest.basis}
\alias{nonest.basis}
\alias{legacy.nonest.basis}
\alias{nonest.basis.qr}
\alias{nonest.basis.matrix}
\alias{nonest.basis.lm}
\alias{nonest.basis.svd}
\alias{nonest.basis.default}
\alias{all.estble}
\alias{is.estble}

\title{Estimability Tools}
\description{
This documents the functions needed to test estimability of linear functions of regression coefficients.
}
\usage{
nonest.basis(x, ...)
\S3method{nonest.basis}{default}(x, ...)
\S3method{nonest.basis}{qr}(x, ...)
\S3method{nonest.basis}{matrix}(x, ...)
\S3method{nonest.basis}{lm}(x, ...)
\S3method{nonest.basis}{svd}(x, tol = 5e-8, ...)

legacy.nonest.basis(x, ...)

all.estble

is.estble(x, nbasis, tol = 1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{For \code{nonest.basis}, an object of a class in \samp{methods("nonest.basis")}. Or, in \code{is.estble}, a numeric vector or matrix for assessing estimability of \samp{sum(x * beta)}, where \code{beta} is the vector of regression coefficients.}
  \item{nbasis}{Matrix whose columns span the null space of the model matrix. Such a matrix is returned by \code{nonest.basis}.}
  \item{tol}{Numeric tolerance for assessing rank or nonestimability. 
  For determining rank, singular values less than \code{tol} times the largest singular value are regarded as zero.
For determining estimability with a nonzero \eqn{x}, \eqn{\beta'x} is assessed by whether or not \eqn{||N'x||^2 < \tau ||x'x||^2}, where \eqn{N} and \eqn{\tau} denote \code{nbasis} and \code{tol}, respectively.}
  \item{\dots}{Additional arguments passed to other methods.}
}

\details{
Consider a linear model \eqn{y = X\beta + E}. If \eqn{X} is not of full rank, it is not possible to estimate \eqn{\beta} uniquely. However, \eqn{X\beta} \emph{is} uniquely estimable, and so is \eqn{a'X\beta} for any conformable vector \eqn{a}. Since \eqn{a'X} comprises a linear combination of the rows of \eqn{X}, it follows that we can estimate any linear function where the coefficients lie in the row space of \eqn{X}. Equivalently, we can check to ensure that the coefficients are orthogonal to the null space of \eqn{X}.

The \code{nonest.basis} method for class \code{'svd'} is not really functional as a method because there is no \code{"svd"} class (at least in R <= 4.2.0). But the function \code{nonest.basis.svd} is exported and may be called directly; it works with results of \code{\link{svd}} or \code{\link{La.svd}}. We \emph{require} \code{x$v} to be the complete matrix of right singular values; but we do not need \code{x$u} at all.

The \code{default} method does serve as an \code{svd} method, in that it only works if \code{x} has the required elements of an SVD result, in which case it passes it to \code{nonest.basis.svd}.
The \code{matrix} method runs \code{nonest.basis.svd(svd(x, nu = 0))}. The \code{lm} method runs the \code{qr} method on \code{x$qr}. 

The function \code{legacy.nonest.basis} is the original default method in early 
versions of the \pkg{estimability} package. It may be called with \code{x} being 
either a matrix or a \code{qr} object, and after obtaining the \code{R} matrix, 
it uses an additional QR decomposition of \code{t(R)} to obtain the needed basis. 
(The current \code{nonest.basis} method for \code{qr} objects is instead based on the 
singular-value decomposition of R, and requires much simpler code.)

The constant \code{all.estble} is simply a 1 x 1 matrix of \code{NA}. This specifies a trivial non-estimability basis, and using it as \code{nbasis} will cause everything to test as estimable. 
}


\value{
When \eqn{X} is not full-rank, the methods for \code{nonest.basis} return a basis for the null space of \eqn{X}. The number of rows is equal to the number of regression coefficients (\emph{including} any \code{NA}s); and the number of columns is equal to the rank deficiency of the model matrix. The columns are orthonormal. If the model is full-rank, then \code{nonest.basis} returns \code{all.estble}. The \code{matrix} method uses \eqn{X} itself, the \code{qr} method uses the \eqn{QR} decomposition of \eqn{X}, and the \code{lm} method recovers the required information from the object.

The function \code{is.estble} returns a logical value (or vector, if \code{x} is a matrix) that is \code{TRUE} if the function is estimable and \code{FALSE} if not. 
}

\references{
Monahan, John F. (2008) \emph{A Primer on Linear Models}, CRC Press. (Chapter 3)
}

\author{
Russell V. Lenth <russell-lenth@uiowa.edu>
}


\examples{
require(estimability)

X <- cbind(rep(1,5), 1:5, 5:1, 2:6)
( nb <- nonest.basis(X) )

SVD <- svd(X, nu = 0)    # we don't need the U part of UDV'
nonest.basis.svd(SVD)    # same result as above

# Test estimability of some linear functions for this X matrix
lfs <- rbind(c(1,4,2,5), c(2,3,9,5), c(1,2,2,1), c(0,1,-1,1))
is.estble(lfs, nb)

# Illustration on 'lm' objects:
warp.lm1 <- lm(breaks ~ wool * tension, data = warpbreaks, 
    subset = -(26:38), 
    contrasts = list(wool = "contr.treatment", tension = "contr.treatment"))
zapsmall(warp.nb1 <- nonest.basis(warp.lm1))

warp.lm2 <- update(warp.lm1, 
    contrasts = list(wool = "contr.sum", tension = "contr.helmert"))
zapsmall(warp.nb2 <- nonest.basis(warp.lm2))

# These bases look different, but they both correctly identify the empty cell
wcells = with(warpbreaks, expand.grid(wool = levels(wool), tension = levels(tension)))
epredict(warp.lm1, newdata = wcells, nbasis = warp.nb1)
epredict(warp.lm2, newdata = wcells, nbasis = warp.nb2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
