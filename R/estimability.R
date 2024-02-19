##############################################################################
#    Copyright (c) 2015-2024 Russell V. Lenth                                #
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

# Obtain an orthonormal basis for nonestimable functions

# Generic
nonest.basis = function(x, ...)
    UseMethod("nonest.basis")


# Legacy code for now-deprecated case of a matrix or qr decomposition
legacy.nonest.basis = function(x, ...) {
    if(!is.qr(x)) {
        if (!is.matrix(x))
            stop("legacy.nonest.basis requires a matrix or qr object")
        x = qr(x)
    }
    rank = x$rank
    tR = t(qr.R(x))
    p = nrow(tR)
    if (rank == p)
        return (all.estble)
    
    # null space of X is same as null space of R in QR decomp
    if (ncol(tR) < p) # add columns if not square
        tR = cbind(tR, matrix(0, nrow=p, ncol=p-ncol(tR)))
    
    # last few rows are zero -- add a diagonal of 1s
    extras = rank + seq_len(p - rank)
    tR[extras, extras] = diag(1, p - rank)
    
    # nbasis is last p - rank cols of Q in QR decomp of tR
    nbasis = qr.Q(qr(tR))[ , extras, drop = FALSE]
    
    # permute the rows via pivot
    nbasis[x$pivot, ] = nbasis
    
    nbasis
}


# Main method -- for class "qr"
# revised to use the svd of R
nonest.basis.qr = function(x, ...) {
    R = qr.R(x)
    cols = which(seq_along(R[1, ]) > x$rank)
    tmp = nbasis = svd(R, nu = 0, nv = ncol(R))$v[, cols, drop = FALSE]
    nbasis[x$pivot, ] = tmp
    nbasis
}

nonest.basis.matrix = function(x, ...)
    nonest.basis.svd(svd(x, nu = 0, nv = ncol(x)), ...)
    ##nonest.basis(qr(x), ...)


nonest.basis.lm = function(x, ...) {
    if (is.null(x$qr))
        x = update(x, method = "qr", qr = TRUE)
    nonest.basis(x$qr)
}

# method for svd class, were it to exist
nonest.basis.svd = function(x, tol = 5e-8, ...) {
    # Note we don't need the 'u' slot at all.
    if(!is.null(x$vt))   # result of La.svd()
        x$v = t(x$vt)
    if (is.null(x$v) || ncol(x$v) < length(x$d))
        stop("We need 'v' to be complete to obtain the basis\n",
             "Run svd() again and exclude the 'nv' argument")
    # we need d to be as long as ncol(v)
    if((deficit <- ncol(x$v) - length(x$d)) > 0)
        x$d = c(x$d, rep(0, deficit))
    w = which(x$d < x$d[1] * tol)
    if (length(w) == 0)
        return(all.estble)
    x$v[, w, drop = FALSE]
}

# default method really designed to suss out an svd() result
nonest.basis.default = function(x, ...) {
    if (!is.null(x$d)) {
        if (!is.null(x$vt) && is.matrix(x$vt))  ## apparently from La.svd
            x$v = t(x$vt)
        if(is.matrix(x$v) && ncol(x$v) == length(x$d))
            return(nonest.basis.svd(x, ...))
    }
    stop("Requires an 'svd()' or 'La.svd()' result")
}



# utility to check estimability of x'beta, given nonest.basis
is.estble = function(x, nbasis, tol = 1e-8) {
    if (is.matrix(x))
        return(apply(x, 1, is.estble, nbasis, tol))
    if(is.na(nbasis[1]))
        TRUE
    else {
        x[is.na(x)] = 0
        chk = as.numeric(crossprod(nbasis, x))
        ssqx = sum(x*x) # BEFORE subsetting x
        # If x really small, don't scale chk'chk
        if (ssqx < tol) ssqx = 1
        sum(chk*chk) < tol * ssqx
    }
}


# nonestimability basis that makes everything estimable
all.estble = matrix(NA)
