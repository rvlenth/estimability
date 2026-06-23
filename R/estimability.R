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

# Obtain an orthonormal basis for nonestimable functions


#' Estimability Tools
#' 
#' This documents the functions needed to test estimability of linear functions 
#' of regression coefficients.
#' 
#' @details Consider a linear model \eqn{y = X\beta + E}. If \eqn{X} is not of full rank,
#' it is not possible to estimate \eqn{\beta} uniquely. However, \eqn{X\beta}
#' \emph{is} uniquely estimable, and so is \eqn{a'X\beta} for any conformable
#' vector \eqn{a}. Since \eqn{a'X} comprises a linear combination of the rows of
#' \eqn{X}, it follows that we can estimate any linear function where the
#' coefficients lie in the row space of \eqn{X}. Equivalently, we can check to
#' ensure that the coefficients are orthogonal to the null space of \eqn{X}.
#' Because estimability constraints are entirely determined by \eqn{X}, these 
#' methods apply to a rich class of models that includes linear models, generalized
#' linear models, mixed models, and models with more general correlated-error
#' structures (e.g., spatial models, time series models).
#'
#' The \code{nonest.basis} method for class \code{'svd'} is not really
#' functional as a method because there is no \code{"svd"} class (at least in R
#' <= 4.2.0). But the function \code{nonest.basis.svd} is exported and may be
#' called directly; it works with results of \code{\link{svd}} or
#' \code{\link{La.svd}}. We \emph{require} \code{x$v} to be the complete matrix
#' of right singular values; but we do not need \code{x$u} at all.
#'
#' The \code{default} method does serve as an \code{svd} method, in that it only
#' works if \code{x} has the required elements of an SVD result, in which case
#' it passes it to \code{nonest.basis.svd}.
#'
#' The \code{matrix} method runs \code{nonest.basis.svd(svd(x, nu = 0))}. The
#' \code{lm} method runs the \code{qr} method on \code{x$qr}.
#'
#' The function \code{legacy.nonest.basis} is the original default method in
#' early versions of the \pkg{estimability} package. It may be called with
#' \code{x} being either a matrix or a \code{qr} object, and after obtaining the
#' \code{R} matrix, it uses an additional QR decomposition of \code{t(R)} to
#' obtain the needed basis. (The current \code{nonest.basis} method for
#' \code{qr} objects is instead based on the singular-value decomposition of R,
#' and requires much simpler code.)
#'
#' The constant \code{all.estble} is simply a 1 x 1 matrix of \code{NA}. This
#' specifies a trivial non-estimability basis, and using it as \code{nbasis}
#' will cause everything to test as estimable.
#'
#' @param x For \code{nonest.basis}, an object of a class in
#'   \samp{methods("nonest.basis")}. Or, in \code{is.estble}, a numeric vector
#'   or matrix for assessing estimability of \samp{sum(x * beta)}, where
#'   \code{beta} is the vector of regression coefficients.
#'
#' @param nbasis Matrix whose columns span the null space of the model matrix.
#'   Such a matrix is returned by \code{nonest.basis}.
#'
#' @param tol Numeric tolerance for assessing rank or nonestimability. For
#'   determining rank, singular values less than \code{tol} times the largest
#'   singular value are regarded as zero. For determining estimability with a
#'   nonzero \eqn{x}, \eqn{\beta'x} is assessed by whether or not \eqn{||N'x||^2
#'   < \tau ||x'x||^2}, where \eqn{N} and \eqn{\tau} denote \code{nbasis} and
#'   \code{tol}, respectively.
#'
#' @param rank Optional integer value. If not \code{NULL}, it is used in place
#'   of \code{tol} to decide the rank of the model matrix
#'
#' @param pivot Optional integer vector of length equal to the number of columns
#'   in the model matrix. If non-\code{NULL}, it is assumed that this specifies
#'   the order of columns in \code{x} relative to the order of predictors in the
#'   model matrix, and the returned non-estimable basis is permuted accordingly.
#'
#' @param ... Additional arguments passed to other methods.
#'
#' @returns When \eqn{X} is not full-rank, the methods for \code{nonest.basis}
#'   return a basis for the null space of \eqn{X}. The number of rows is equal
#'   to the number of regression coefficients (\emph{including} any \code{NA}s);
#'   and the number of columns is equal to the rank deficiency of the model
#'   matrix. The columns are orthonormal. If the model is full-rank, then
#'   \code{nonest.basis} returns \code{all.estble}. The \code{matrix} method
#'   uses \eqn{X} itself, the \code{qr} method uses the \eqn{QR} decomposition
#'   of \eqn{X}, and the \code{lm} method recovers the required information from
#'   the object.
#'
#'   The function \code{is.estble} returns a logical value (or vector, if
#'   \code{x} is a matrix) that is \code{TRUE} if the function is estimable and
#'   \code{FALSE} if not.
#'
#' @section Use with Cholesky decompositions: If the \code{\link{chol}} function
#'   is called with \code{pivot = TRUE}, the returned matrix is suitable to
#'   provide to the \code{matrix} method of \code{nonest.basis}, and the
#'   returned basis is correctly pivoted without specifying \code{pivot}
#'   explicitly.
#'
#'   This works because the Cholesky decomposition \eqn{R} of \eqn{X'X} (as constructed by
#'   \code{chol()}) satisfies \eqn{X'X = R'R}, implying that the rows of
#'   \eqn{X'X} are linear combinations of the rows of \eqn{R}; and they are both
#'   of equal rank - so \eqn{R} and \eqn{X'X} have the same row space. 
#'
#' @order 1
#' @export
#'
#' @examples
#' require(estimability)
#'
#' X <- cbind(rep(1,5), 1:5, 5:1, 2:6)
#' ( nb <- nonest.basis(X) )
#'
#' # via the singular-value decomposition
#' SVD <- svd(X, nu = 0)    # we don't need the U part of UDV'
#' nonest.basis.svd(SVD)    # same result as above
#'
#' # Via the Cholesky decomposition of X'X
#' ch <- chol(t(X) %*% X, pivot = TRUE) |> suppressWarnings()
#' ( nb.chol <- nonest.basis(ch) )
#'
#' # Test estimability of some linear functions for this X matrix
#' lfs <- rbind(c(1,4,2,5), c(2,3,9,5), c(1,2,2,1), c(0,1,-1,1))
#' is.estble(lfs, nb)
#' is.estble(lfs, nb.chol) # same results even though nb.chol is not identical to nb
#'
#' # Illustration on 'lm' objects:
#' warp.lm1 <- lm(breaks ~ wool * tension, data = warpbreaks,
#'     subset = -(26:38),
#'     contrasts = list(wool = "contr.treatment", tension = "contr.treatment"))
#' zapsmall(warp.nb1 <- nonest.basis(warp.lm1))
#'
#' # different parameterization of the same model:
#' warp.lm2 <- update(warp.lm1,
#'     contrasts = list(wool = "contr.sum", tension = "contr.helmert"))
#' zapsmall(warp.nb2 <- nonest.basis(warp.lm2))
#'
#' # These bases depend of the parameterizations, but they both 
#' # correctly identify the empty cell
#' wcells = with(warpbreaks, expand.grid(wool = levels(wool), tension = levels(tension)))
#' epredict(warp.lm1, newdata = wcells, nbasis = warp.nb1)
#' epredict(warp.lm2, newdata = wcells, nbasis = warp.nb2)
#' 
nonest.basis = function(x, ...)
    UseMethod("nonest.basis")


# Legacy code for now-deprecated case of a matrix or qr decomposition
#' @rdname nonest.basis
#' @order 25
#' @export
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
#' @rdname nonest.basis
#' @order 11
#' @export
nonest.basis.qr = function(x, ...) {
    R = qr.R(x)
    cols = which(seq_along(R[1, ]) > x$rank)
    tmp = nbasis = svd(R, nu = 0, nv = ncol(R))$v[, cols, drop = FALSE]
    nbasis[x$pivot, ] = tmp
    nbasis
}

#' @rdname nonest.basis
#' @order 12
#' @export
nonest.basis.matrix = function(x, tol = 5e-8, rank = attr(x, "rank"), pivot = attr(x, "pivot"), ...) {
    nonest.basis.svd(svd(x, nu = 0, nv = ncol(x)), rank = rank, pivot = pivot, ...)
}

#' @rdname nonest.basis
#' @order 13
#' @export
nonest.basis.lm = function(x, ...) {
    if (is.null(x$qr))
        x = update(x, method = "qr", qr = TRUE)
    nonest.basis(x$qr)
}

# method for svd class, were it to exist
# If non-null, rank is used in lieu of tol, and results are permuted using pivot
#' @rdname nonest.basis
#' @order 14
#' @export nonest.basis.svd
nonest.basis.svd = function(x, tol = 5e-8, rank = NULL, pivot = NULL, ...) {
    # Note we don't need the 'u' slot at all.
    if(!is.null(x$vt))   # result of La.svd()
        x$v = t(x$vt)
    if (is.null(x$v) || ncol(x$v) < length(x$d))
        stop("We need 'v' to be complete to obtain the basis\n",
             "Run svd() again and exclude the 'nv' argument")
    if(is.null(rank)) {
        if((deficit <- ncol(x$v) - length(x$d)) > 0) # we need d to be as long as ncol(v)
            x$d = c(x$d, rep(0, deficit))
        w = which(x$d < x$d[1] * tol)
    }
    else
        w = seq_len(ncol(x$v) - rank) + rank
    if (length(w) == 0)
        return(all.estble)
    bas = x$v[, w, drop = FALSE]
    if(!is.null(pivot))
        bas = bas[order(pivot), , drop = FALSE]
    bas
}

# default method really designed to suss out an svd() result
#' @rdname nonest.basis
#' @order 10
#' @export
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
#' @rdname nonest.basis
#' @order 51
#' @export
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
# Note without the @format, roxygen2 thinks this is a data set
#' @rdname nonest.basis
#' @format NULL
#' @order 52
#' @export
all.estble = matrix(NA)
