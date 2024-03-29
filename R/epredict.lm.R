##############################################################################
#    Copyright (c) 2015-2016 Russell V. Lenth                                #
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
epredict = function(object, ...)
    UseMethod("epredict")

epredict.lm = function(object, newdata, ..., 
        type = c("response", "terms", "matrix", "estimability"), 
        nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)

epredict.glm = function(object, newdata, ..., 
        type = c("link", "response", "terms", "matrix", "estimability"), 
        nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)

epredict.mlm = function(object, newdata, ..., 
            type = c("response", "matrix", "estimability"), 
            nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, type[1], nonest.tol, nbasis, ...)


# Generic for eupdate -- adds nonest basis to object
eupdate = function(object, ...)
    UseMethod("eupdate")

eupdate.lm = function(object, ...) {
    if (length(list(...)) > 0)
        object = do.call("update", list(object = object, ...))
    object$nonest = nonest.basis(object)
    object
}
