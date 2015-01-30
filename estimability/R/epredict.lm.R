# Patch for predict.lm, predict.glm, predict.mlm
#
# If newdata present, and fit is rank-deficient, 
# we check estimability and replace any non-est cases with NA
# Use options(estimability.quiet = TRUE) to suppress message
# Use options(estimability.suppress = TRUE) to override this patch

# Convenience utility - return 
#    TRUE  if an option is present and TRUE, 
#    FALSE if option is NULL or FALSE
.falseIfNull = function(option)
    !is.null(opt <- getOption(option)) && opt


# Main workhorse -- call with stats-library predict function
.patch.predict = function(object, newdata, 
        nonest.tol = 1e-8, nbasis = object$nonest, ...) {
    
    if(missing(newdata)) 
        predict(object = object, ...)
    
    else {
        if (all(!is.na(object$coefficients)))
            return (predict(object = object, newdata = newdata, ...))
        
        w.handler <- function(w){ # suppress the incorrect warning
            if (!is.na(pmatch("prediction from a rank-deficient", w$message)))
                invokeRestart("muffleWarning")
        }
        result = withCallingHandlers(predict(object = object, newdata = newdata, ...),
            warning = w.handler)
        
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
        if (any(nonest)) {
            if (is.matrix(result))
                result[nonest, ] = NA
            else if (is.list(result)) {
                result$fit[nonest] = NA
                result$se.fit[nonest] = NA
            }
            else
                result[nonest] = NA
            if(!.falseIfNull("estimability.quiet"))
                message("Non-estimable cases are replaced by 'NA'")
        }
        
        result
    }
}


# Generic for epredict
epredict = function(object, ...)
    UseMethod("epredict")

epredict.lm = function(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, nonest.tol, nbasis, ...)

epredict.glm = function(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, nonest.tol, nbasis, ...)

epredict.mlm = function(object, newdata, ..., nonest.tol = 1e-8, nbasis = object$nonest)
    .patch.predict(object, newdata, nonest.tol, nbasis, ...)
