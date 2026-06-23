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

# Obtain an estimable subspace from the rows of a matrix L
# i.e., B %*% L such that B %*% L %*% N = 0 (where N = nbasis)
# Thus, (LN)'B' = 0, i.e., B' is in null space of LN'
# We are tooled-up to find that!
#
# The function returns BL, with B as an attribute



#' Find an estimable subspace
#' 
#' Determine a transformation \code{B} of the rows of a matrix \code{L}
#' such that \code{B \%*\% L} is estimable.
#' A practical example is in jointly testing a set of contrasts \code{L}
#' in a linear model, and we need to restrict to the subspace spanned by
#' the rows of \code{L} that are estimable.
#' 
#' @details We require \code{B} such that all the rows of \code{M = B \%*\% L}
#' are estimable, i.e. orthogonal to the columns of \code{nbasis}.
#' Thus, we need \code{B \%*\% L \%*\% nbasis} to be zero, or equivalently,
#' \code{t(B)} must be in the null space of \code{t(L \%*\% nbasis)}.
#' This can be found using \code{\link{nonest.basis}}.
#'
#' @param L A matrix of dimensions \emph{k} by \emph{p}
#' @param nbasis A \emph{k} by \emph{b} matrix whose columns form a
#'   basis for non-estimable linear functions -- such as is returned 
#'   by \code{\link{nonest.basis}}
#' @param tol Numeric tolerance for assessing nonestimability. See
#'   \code{\link{is.estble}}.
#'
#' @returns
#'   An \emph{r} by \emph{p} matrix \code{M = B \%*\% L} 
#'   whose rows are all orthogonal to the columns of 
#'   \code{nbasis}. The matrix \code{B} is attached as \code{attr(M, "B")}.
#'   Note that if any rows of \code{L} were non-estimable, then \emph{r}
#'   will be less than \emph{k}. In fact, if there are no estimable functions 
#'   in the row space of \code{L}, then \emph{r} = 0.
#' 
#' @export
#'
#' @examples
#' ### Find a set of estimable interaction contrasts for a 3 x 4 design 
#' ### with two empty cells.
#' des <- expand.grid(A = factor(1:3), B = factor(1:4))
#' des <- des[-c(5, 12), ]  # cells (2,2) and (3,4) are empty
#' 
#' X <- model.matrix(~ A * B, data = des)
#' N <- nonest.basis(X)
#' 
#' L <- cbind(matrix(0, nrow = 6, ncol = 6), diag(6))
#' # i.e., give nonzero weight only to interaction effects
#' 
#' estble.subspace(L, N)
#' 
#' # Tougher demo: create a variation where all rows of L are non-estimable
#' set.seed(1921.0923)  # Franklin Graybill's birthday
#' LL <- matrix(rnorm(36), ncol = 6) %*% L
#' estble.subspace(LL, N)
#' 
estble.subspace = function(L, nbasis, tol = 1e-8) {
    if (all(apply(L, 1, is.estble, nbasis, tol)))
        B = diag(nrow(L))
    else {
        LN = L %*% nbasis
        LN[abs(LN) <= tol] = 0   # don't be jerked around by small values
        B = t(nonest.basis(t(LN)))
    }
    if (is.na(B[1])) # nothing is estimable
        result = matrix(0, nrow = 0, ncol = ncol(L))
    else
        result = B %*% L
    attr(result, "B") = B
    result
}