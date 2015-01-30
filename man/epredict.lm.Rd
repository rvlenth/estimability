\name{epredict}
\alias{epredict}
\alias{epredict.lm}
\alias{epredict.glm}
\alias{epredict.mlm}

\title{
Predictions with estimability checking
}
\description{
These functions call the corresponding S3 \code{predict} methods in the \pkg{stats} package, but with a check for estimability of new predictions, and with appropriate actions for non-estimable cases.
}
\usage{
\S3method{epredict}{lm}(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
\S3method{epredict}{glm}(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
\S3method{epredict}{mlm}(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
}
\arguments{
  \item{object}{An object ingeriting from \code{lm}}
  \item{newdata}{A \code{data.frame} containing predictor combinations for new predictions}
  \item{\dots}{Arguments passed to \code{\link{predict}}}
  \item{nonest.tol}{Tolerance used by \code{\link{is.estble}} to check estimability of new predictions}
  \item{nbasis}{Basis for the row-null space, e.g., a result of a call to \code{\link{nonest.basis}}. If \code{nbasis} is \code{NULL}, a basis is constructed from \code{object}.}
}
\details{
If \code{newdata} is missing or \code{object} is not rank-deficient, this method passes its arguments directly to the same method in the \pkg{stats} library. In rank-deficient cases with \code{newdata} provided, each row of \code{newdata} is tested for estimability against the null basis provided in \code{nbasis}. Any non-estimable cases found are replaced with \code{NA}s.

If you anticipate making several \code{epredict} calls with new data, it improves efficiency to either obtain the null basis and provide it in the call, or add it to \code{object} with the name \code{"nonest"}. See the example below.
}

\value{
The same as the result of a call to the \code{predict} method in the \pkg{stats} package, except rows or elements corresponding to non-estimable predictor combinations are set to \code{NA}.}

\author{
Russell V. Lenth <russell-lenth@uiowa.edu>
}

\note{
The usual rank-deficiency warning from \code{stats::predict} is suppressed; but when non-estimable cases are found, a message is displayed explaining that these results were replaced by \code{NA}. If you wish that message suppressed, use \samp{options(estimability.quiet = TRUE)}.
}


\seealso{
  \code{\link[stats]{predict.lm}} in the \pkg{stats} package;   
  \code{\link{nonest.basis}}.
}
\examples{
require("estimability")

# Fake data where x3 and x4 depend on x1, x2, and intercept
x1 <- -4:4
x2 <- c(-2,1,-1,2,0,2,-1,1,-2)
x3 <- 3*x1 - 2*x2
x4 <- x2 - x1 + 4
y <- 1 + x1 + x2 + x3 + x4 + c(-.5,.5,.5,-.5,0,.5,-.5,-.5,.5)

# Different orderings of predictors produce different solutions
mod1234 <- lm(y ~ x1 + x2 + x3 + x4)
mod4321 <- lm(y ~ x4 + x3 + x2 + x1)
 
# test data:
testset <- data.frame( 
              x1 = c(3, 6, 6, 0, 0, 1), 
              x2 = c(1, 2, 2, 0, 0, 2), 
              x3 = c(7,14,14, 0, 0, 3), 
              x4 = c(2, 4, 0, 4, 0, 4))

# Look at predictions when we don't check estimability
suppressWarnings( # Disable the warning from stats::predict.lm
    rbind(p1234 = predict(mod1234, newdata = testset),
          p4321 = predict(mod4321, newdata = testset)))

# Compare with results when we do check:
rbind(p1234 = epredict(mod1234, newdata = testset),
      p4321 = epredict(mod4321, newdata = testset))

# Note that estimable cases have the same predictions

# To prevent re-constructing the basis in future calls:
mod1234$nonest <- nonest.basis(mod1234$qr)

# When row spaces are the same, bases are interchangeable
# so long as you account for the ordering of parameters:
epredict(mod4321, newdata = testset, 
    nbasis = mod1234$nonest[c(1,5:2), ])

}

\keyword{ models }
\keyword{ regression }
