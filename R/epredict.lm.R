##############################################################################
#    Copyright (c) 2015-2026 Russell V. Lenth                                #
#                                                                            #
#    This file is part of the estimability package for R (*estimability*)    #
#                                                                            #
#    *estimability* is free software: you can redistribute it and/or modify  #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *estimability* is distributed in the hope that it will be useful,       #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    A copy of the GNU General Public License is available at                #
#    <https://www.r-project.org/Licenses/>                                   #
##############################################################################

# Patch for predict.lm, predict.glm, predict.mlm
#
# If newdata present, and fit is rank-deficient, 
# we check estimability and replace any non-est cases with NA
# Use options(estimability.quiet = TRUE) to suppress message
# Use options(estimability.suppress = TRUE) to override this patch

# Main workhorse -- call with stats-library predict function
.patch.predict = function(object, newdata, type,
        nonest.tol = 1e-8, nbasis = object$nonest, ...) {
    
    if(missing(newdata)) 
        predict(object = object, type = type, ...)
    
    else {
        type = match.arg(type, c("response", "link", "terms", "matrix", "estimability"))
        if (all(!is.na(object$coefficients)) && (type != "matrix"))
            if (type == "estimability")
                return (rep(TRUE, nrow(newdata)))
            else
                return (predict(object = object, newdata = newdata, type = type, ...))
        
        if(is.null(nbasis)) {
            if (!is.null(qr <- object$qr))
                nbasis = nonest.basis(qr)
            else
                nbasis = nonest.basis(model.matrix(object))
        }
        trms = delete.response(terms(object))
        m = model.frame(trms, newdata, na.action = na.pass, xlev = object$xlevels)
        X = model.matrix(trms, m, contrasts.arg = object$contrasts)
        nonest = !is.estble(X, nbasis, nonest.tol)
        
        if (type == "estimability")
            return (!nonest)
        else if (type == "matrix") {
            attr(X, "estble") = !nonest
            return (X)
        }
        # (else) we have a type anticipated by stats::predict
        
        w.handler <- function(w){ # suppress the incorrect warning
            if (!is.na(pmatch("prediction from a rank-deficient", w$message)))
                invokeRestart("muffleWarning")
        }
        result = withCallingHandlers(
            suppressWarnings(predict(object = object, newdata = newdata, type = type, 
                    rankdeficient = "simple", ...)),
            warning = w.handler)
        
        if (any(nonest)) {
            if (is.matrix(result))
                result[nonest, ] = NA
            else if (is.list(result)) {
                result$fit[nonest] = NA
                result$se.fit[nonest] = NA
            }
            else
                result[nonest] = NA
            if(getOption("estimability.verbose", FALSE))
                message("Note: Non-estimable cases were replaced by 'NA'")
        }
        
        result
    }
}


# Generic for epredict
#' Estimability Enhancements for \code{lm} and Relatives
#' 
#' These functions call the corresponding S3 \code{predict} methods in the 
#' \pkg{stats} package, but with a check for estimability of new predictions, 
#' and with appropriate actions for non-estimable cases.
#'
#' If \code{newdata} is missing or \code{object} is not rank-deficient, this method 
#' irectly to the same method in the \pkg{stats} library. In rank-deficient cases with 
#' \code{newdata} provided, each row of \code{newdata} is tested for estimability against 
#' the null basis provided in \code{nbasis}. Any non-estimable cases found are replaced with \code{NA}s.
#' 
#' The \code{type} argument is passed to \code{\link[stats]{predict}} 
#' when it is one of \code{"response"}, \code{"link"}, or \code{"terms"}. 
#' With \code{newdata} present and \code{type = "matrix"}, the model 
#' matrix for \code{newdata} is returned, with an attribute \code{"estble"} 
#' that is a logical vector of length \samp{nrow(newdata)} indicating whether 
#' each row is estimable. With \code{type = "estimability"}, just the logical 
#' vector is returned.
#' 
#' If you anticipate making several \code{epredict} calls with new data, 
#' it improves efficiency to either obtain the null basis and provide it 
#' in the call, or add it to \code{object} with the name \code{"nonest"} 
#' (perhaps via a call to \code{eupdate}).
#' 
#' \code{eupdate} is an S3 generic function with a method provided for \code{"lm"} 
#' It updates the object according to any arguments in \code{...}, then obtains 
#' the updated object's nonestimable basis and returns it in \code{object$nonest}.

#' @param object An object inheriting from \code{lm}
#' @param ... Arguments passed to \code{\link{predict}} or \code{\link{update}}
#' @param newdata A \code{data.frame} containing predictor combinations for new predictions
#' @param nonest.tol Tolerance used by \code{\link{is.estble}} to check estimability of new predictions
#' @param type haracter string specifying the desired result. See Details.
#' @param nbasis a basis for the null space, e.g., a result of a call to \code{\link{nonest.basis}}. If \code{nbasis} is \code{NULL}, a basis is constructed from \code{object}.
#'
#' @returns The same as the result of a call to the \code{predict} method in the 
#'   \pkg{stats} package, except rows or elements corresponding to non-estimable 
#'   predictor combinations are set to \code{NA}. 
#'   The value for \code{type} is  \code{"matrix"} or \code{"estimability"} 
#'   is explained under details.
#'   
#' @note
#' The capabilities of the \code{epredict} function for \code{lm} objects is provided 
#' by \code{\link[stats]{predict.lm}} (if using R version 4.3.0 or later) with \code{rankdeficient = "NA"};
#' however, \code{epredict} 
#' uses \pkg{estimability}'s own criteria to determine which predictions
#' are set to \code{NA}. An advantage of using \code{epredict} is one of
#' efficiency: we can compute the null basis once and for all and have it available
#' additional predictions, whereas \code{predict.lm} will re-compute it each time.
#' If the user wishes to see a message explaining why \code{NA}s were displayed,
#' set \samp{options(estimability.verbose = TRUE)}.
#' 
#' @order 1 
#' @export
#' 
#' @importFrom stats delete.response
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats na.pass
#' @importFrom stats predict
#' @importFrom stats terms
#' @importFrom stats update
#'
#' @examples
#' require("estimability")
#' 
#' # Fake data where x3 and x4 depend on x1, x2, and intercept
#' x1 <- -4:4
#' x2 <- c(-2,1,-1,2,0,2,-1,1,-2)
#' x3 <- 3*x1 - 2*x2
#' x4 <- x2 - x1 + 4
#' y <- 1 + x1 + x2 + x3 + x4 + c(-.5,.5,.5,-.5,0,.5,-.5,-.5,.5)
#' 
#' # Different orderings of predictors produce different solutions
#' mod1234 <- lm(y ~ x1 + x2 + x3 + x4)
#' mod4321 <- eupdate(lm(y ~ x4 + x3 + x2 + x1))
#' # (Estimability checking with mod4321 will be more efficient because
#' #  it will not need to recreate the basis)
#' mod4321$nonest
#' 
#'  
#' # test data:
#' testset <- data.frame( 
#'               x1 = c(3,  6,  6,  0,  0,  1), 
#'               x2 = c(1,  2,  2,  0,  0,  2), 
#'               x3 = c(7, 14, 14,  0,  0,  3), 
#'               x4 = c(2,  4,  0,  4,  0,  4))
#' 
#' # Look at predictions when we don't check estimability
#' suppressWarnings( # Disable the warning from stats::predict.lm
#'     rbind(p1234 = predict(mod1234, newdata = testset),
#'           p4321 = predict(mod4321, newdata = testset)))
#' 
#' # Compare with results when we do check:
#' rbind(p1234 = epredict(mod1234, newdata = testset),
#'       p4321 = epredict(mod4321, newdata = testset))
#' 
#' # stats::predict has same capability for lm objects starting in version 4.3.0:
#' if((R.Version()$major >= 4) && (R.Version()$minor >= 3))
#'   stats::predict(mod1234, newdata = testset, rankdeficient = "NA")
#' 
#' # Note that estimable cases have the same predictions
#' 
#' # change mod1234 and include nonest basis 
#' mod134 <- eupdate(mod1234, . ~ . - x2, subset = -c(3, 7))
#' mod134$nonest
#' 
#' # When row spaces are the same, bases are interchangeable
#' # so long as you account for the ordering of parameters:
#' epredict(mod4321, newdata = testset, type = "estimability",
#'     nbasis = nonest.basis(mod1234)[c(1,5:2), ])
#'     
#' # Comparison with predict.lm in R >= 4.3.0
#' if((R.Version()$major >= 4) && (R.Version()$minor >= 3))
#'   stats::predict(mod4321, newdata = testset, rankdeficient = "NA")
#' 
#' \dontrun{
#' ### Additional illustration
#' example(nonest.basis)  ## creates model objects warp.lm1 and warp.lm2
#' 
#' # The two models have different contrast specs. But the empty cell
#' # is correctly identified in both:
#' fac.cmb <- expand.grid(wool = c("A", "B"), tension = c("L", "M", "H"))
#' cbind(fac.cmb, 
#'       pred1 = epredict(warp.lm1, newdata = fac.cmb), 
#'       pred2 = epredict(warp.lm2, newdata = fac.cmb))
#' }
#' 
epredict = function(object, ...) {
    UseMethod("epredict")
}

#' @rdname epredict
#' @order 2
#' @export
epredict.lm = function(object, newdata, ..., 
        type = c("response", "terms", "matrix", "estimability"), 
        nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)

#' @rdname epredict
#' @order 3
#' @export
epredict.glm = function(object, newdata, ..., 
        type = c("link", "response", "terms", "matrix", "estimability"), 
        nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)

#' @rdname epredict
#' @order 4
#' @export
epredict.mlm = function(object, newdata, ..., 
            type = c("response", "matrix", "estimability"), 
            nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)


# Generic for eupdate -- adds nonest basis to object
#' @rdname epredict
#' @order 21
#' @export
eupdate = function(object, ...)
    UseMethod("eupdate")

#' @rdname epredict
#' @order 22
#' @export
eupdate.lm = function(object, ...) {
    if (length(list(...)) > 0)
        object = do.call("update", list(object = object, ...))
    object$nonest = nonest.basis(object)
    object
}
