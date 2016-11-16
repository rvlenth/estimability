%    Copyright (c) 2015-2016 Russell V. Lenth                                #

\name{epredict}
\alias{epredict}
\alias{epredict.lm}
\alias{epredict.glm}
\alias{epredict.mlm}
\alias{eupdate}
\alias{eupdate.lm}

\title{
Estimability Enhancements for \code{lm} and Relatives
}
\description{
These functions call the corresponding S3 \code{predict} methods in the \pkg{stats} package, but with a check for estimability of new predictions, and with appropriate actions for non-estimable cases.
}
\usage{
\S3method{epredict}{lm}(object, newdata, ..., 
    type = c("response", "terms", "matrix", "estimability"), 
    nonest.tol = 1e-8, nbasis = object$nonest)
\S3method{epredict}{glm}(object, newdata, ..., 
    type = c("link", "response", "terms", "matrix", "estimability"), 
    nonest.tol = 1e-8, nbasis = object$nonest)
\S3method{epredict}{mlm}(object, newdata, ..., 
    type = c("response", "matrix", "estimability"), 
    nonest.tol = 1e-8, nbasis = object$nonest)
    
eupdate(object, ...)    
}
\arguments{
  \item{object}{An object ingeriting from \code{lm}}
  \item{newdata}{A \code{data.frame} containing predictor combinations for new predictions}
  \item{\dots}{Arguments passed to \code{\link{predict}} or \code{\link{update}}}
  \item{nonest.tol}{Tolerance used by \code{\link{is.estble}} to check estimability of new predictions}
  \item{type}{Character string specifying the desired result. See Details.}
  \item{nbasis}{Basis for the null space, e.g., a result of a call to \code{\link{nonest.basis}}. If \code{nbasis} is \code{NULL}, a basis is constructed from \code{object}.}
}
\details{
If \code{newdata} is missing or \code{object} is not rank-deficient, this method passes its arguments directly to the same method in the \pkg{stats} library. In rank-deficient cases with \code{newdata} provided, each row of \code{newdata} is tested for estimability against the null basis provided in \code{nbasis}. Any non-estimable cases found are replaced with \code{NA}s.

The \code{type} argument is passed to \code{\link[stats]{predict}} when it is one of \code{"response"}, \code{"link"}, or \code{"terms"}. With \code{newdata} present and \code{type = "matrix"}, the model matrix for \code{newdata} is returned, with an attribute \code{"estble"} that is a logical vector of length \samp{nrow(newdata)} indicating whether each row is estimable. With \code{type = "estimability"}, just the logical vector is returned.

If you anticipate making several \code{epredict} calls with new data, it improves efficiency to either obtain the null basis and provide it in the call, or add it to \code{object} with the name \code{"nonest"} (perhaps via a call to \code{eupdate}).

\code{eupdate} is an S3 generic function with a method provided for \code{"lm"} objects. It updates the object according to any arguments in \code{...}, then obtains the updated object's nonestimable basis and returns it in \code{object$nonest}.
}

\value{
The same as the result of a call to the \code{predict} method in the \pkg{stats} package, except rows or elements corresponding to non-estimable predictor combinations are set to \code{NA}. The value for \code{type} is  \code{"matrix"} or \code{"estimability"} is explained under details.}

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
mod4321 <- eupdate(lm(y ~ x4 + x3 + x2 + x1))
# (Estimability checking with mod4321 will be more efficient because
#  it will not need to recreate the basis)
mod4321$nonest

 
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

# Compare with results when we do check:
rbind(p1234 = epredict(mod1234, newdata = testset),
      p4321 = epredict(mod4321, newdata = testset))

# Note that estimable cases have the same predictions

# change mod1234 and include nonest basis 
mod134 <- eupdate(mod1234, . ~ . - x2, subset = -c(3, 7))
mod134$nonest

# When row spaces are the same, bases are interchangeable
# so long as you account for the ordering of parameters:
epredict(mod4321, newdata = testset, type = "estimability",
    nbasis = nonest.basis(mod1234)[c(1,5:2), ])
    
\dontrun{
### Additional illustration
example(nonest.basis)  ## creates model objects warp.lm1 and warp.lm2

# The two models have different contrast specs. But the empty cell
# is correctly identified in both:
fac.cmb <- expand.grid(wool = c("A", "B"), tension = c("L", "M", "H"))
cbind(fac.cmb, 
      pred1 = epredict(warp.lm1, newdata = fac.cmb), 
      pred2 = epredict(warp.lm2, newdata = fac.cmb))
} % end of \dontrun   

}

\keyword{ models }
\keyword{ regression }
