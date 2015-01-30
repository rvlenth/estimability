\name{nonest.basis}
\alias{nonest.basis}
\alias{all.estble}
\alias{is.estble}

\title{Estimability tools}
\description{
This documents the functions needed to test estimability of linear functions of regression coefficients.
}
\usage{
nonest.basis(qrX)
all.estble

is.estble(x, nbasis, tol = 1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qrX}{A model matrix, or the QR decomposition thereof (from \code{\link{qr}} with \samp{LAPACK=FALSE}). The latter is preferred if already available, as it saves computation.
  }
  \item{x}{Numeric vector or matrix for assessing estimability of \samp{sum(x * beta)}, where \code{beta} is the vector of regression coefficients.}
  \item{nbasis}{Matrix whose columns span the null-row space of the model matrix. Such a matrix is returned by \code{nonest.basis}.}
  \item{tol}{Numeric tolerance for assessing noestimability. For nonzero \eqn{x}, estimability of \eqn{\beta'x} is assessed by whether or not
  \eqn{||N'x||^2 < \tau ||x'x||^2}, where \eqn{N} and \eqn{\tau} denote \code{nbasis} and \code{tol}, respectively.}
}

\details{
Consider a linear model \eqn{y = X\beta + E}. If \eqn{X} is not of full rank, it is not possible to estimate \eqn{\beta} uniquely. However, \eqn{X\beta} \emph{is} uniquely estimable, and so is \eqn{a'X\beta} for any conformable vector \eqn{a}. Since \eqn{a'X} comprises a linear combination of the rows of \eqn{X}, it follows that we can estimate any linear function where the coefficients lie in the row space of \eqn{X}. Equivalently, we can check to ensure that the coefficients are orthogonal to the null space of \eqn{X'}. 

The constant \code{all.estble} is simply a 1 x 1 matrix of \code{NA}. This specifies a \sQuote{null} non-estimability basis, and using it as \code{nbasis} will cause everything to test as estimable. 
}


\value{
When \eqn{X} is not full-rank, the function \code{nonest.basis} returns a basis for the null space of \eqn{X'}. The number of rows is equal to \samp{length(coef(object))}; and the number of columns is equal to the rank deficiency of the model matrix. The columns are orthonormal. If the model is full-rank, then \code{nonest.basis} returns \code{all.estble}.

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

lfs <- rbind(c(1,4,2,5), c(2,3,9,5), c(1,2,2,1), c(0,1,-1,1))
is.estble(lfs, nb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
