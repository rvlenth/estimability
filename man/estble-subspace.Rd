%    Copyright (c) 2015-2018 Russell V. Lenth

\name{estble.subspace}
\alias{estble.subspace}

\title{Find an estimable subspace}

\description{
  Determine a transformation \code{B} of the rows of a matrix \code{L}
  such that \code{B \%*\% L} is estimable.
  A practical example is in jointly testing a set of contrasts \code{L}
  in a linear model, and we need to restrict to the subspace spanned by
  the rows of \code{L} that are estimable.
}

\usage{
estble.subspace (L, nbasis, tol = 1e-8)
}

\arguments{
  \item{L}{A matrix of dimensions \emph{k} by \emph{p}}
  \item{nbasis}{A \emph{k} by \emph{b} matrix whose columns form a
    basis for non-estimable linear functions -- such as is returned 
    by \code{\link{nonest.basis}}}
  \item{tol}{Numeric tolerance for assessing nonestimability. See
    \code{\link{is.estble}}.}
}

\details{
  We require \code{B} such that all the rows of \code{M = B \%*\% L}
  are estimable, i.e. orthogonal to the columns of \code{nbasis}.
  Thus, we need \code{B \%*\% L \%*\% nbasis} to be zero, or equivalently,
  \code{t(B)} must be in the null space of \code{t(L \%*\% nbasis)}.
  This can be found using \code{\link{nonest.basis}}.
}


\value{
  An \emph{r} by \emph{p} matrix \code{M = B \%*\% L} 
  whose rows are all orthogonal to the columns of 
  \code{nbasis}. The matrix \code{B} is attached as \code{attr(M, "B")}.
  Note that if any rows of \code{L} were non-estimable, then \emph{r}
  will be less than \emph{k}. In fact, if there are no estimable functions 
  in the row space of \code{L}, then \emph{r} = 0.
}


\author{
Russell V. Lenth <russell-lenth@uiowa.edu>
}


\examples{
### Find a set of estimable interaction contrasts for a 3 x 4 design 
### with two empty cells.
des <- expand.grid(A = factor(1:3), B = factor(1:4))
des <- des[-c(5, 12), ]  # cells (2,2) and (3,4) are empty

X <- model.matrix(~ A * B, data = des)
N <- nonest.basis(X)

L <- cbind(matrix(0, nrow = 6, ncol = 6), diag(6))
# i.e., give nonzero weight only to interaction effects

estble.subspace(L, N)

# Tougher demo: create a variation where all rows of L are non-estimable
LL <- matrix(rnorm(36), ncol = 6) \%*\% L
estble.subspace(LL, N)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
