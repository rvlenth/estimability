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

# Obtain an orthonormal basis for nonestimable functions

# Generic
nonest.basis = function(x, ...)
    UseMethod("nonest.basis")


# Main method -- for class "qr"
nonest.basis.qr = function(x, ...) {
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


nonest.basis.matrix = function(x, ...)
    nonest.basis(qr(x), ...)


nonest.basis.lm = function(x, ...) {
    if (is.null(x$qr))
        x = update(x, method = "qr", qr = TRUE)
    nonest.basis(x$qr)
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
