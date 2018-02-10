##############################################################################
#    Copyright (c) 2015-2018 Russell V. Lenth                                #
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
estble.subspace = function(L, nbasis, tol = 1e-8) {
    if (all(apply(L, 1, is.estble, nbasis, tol)))
        B = diag(nrow(L))
    else {
        LN = L %*% nbasis
        B = t(nonest.basis(t(LN)))
    }
    if (is.na(B[1])) # nothing is estimable
        result = L[-seq_len(nrow(L)), , drop = FALSE]  # nrow = 0
    else
        result = B %*% L
    attr(result, "B") = B
    result
}