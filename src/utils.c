#include <Rinternals.h>
#include <htslib/kfunc.h>

/*! @function
 @abstract Simple wrapper to run fisher exact test using htslib
 @param  mat   integer matrix, with 4 rows per column,
               populated with base counts representing a_ref, a_alt, b_ref, b_alt

 @return
  Numeric vector with p-values
 */
SEXP fisher_exact(SEXP mat) {
    int nr = 4;
    if (!Rf_isMatrix(mat) || !(Rf_nrows(mat) == nr)) {
        Rf_error("'mat' must be matrix with 4 rows");
    }

    if (!Rf_isInteger(mat)) {
       Rf_error("'mat' must be an integer matrix");
    }


    int nc = Rf_ncols(mat);
    SEXP res = PROTECT(Rf_allocVector(REALSXP, nc));
    int n11, n12, n21, n22;

    // R matrices are stored as vectors in column major order
    // approach for iteration based on base R colSums C source (do_colsum in array.c )
    for (int j = 0; j < nc; j++) {
        int *ix = INTEGER(mat) + (R_xlen_t)nr*j;
        n11 = *ix++;
        n12 = *ix++;
        n21 = *ix++;
        n22 = *ix;
        double left, right, fisher;
        kt_fisher_exact(n11,n12,n21,n22, &left,&right,&fisher);
        REAL(res)[j] = fisher;
    }

    UNPROTECT(1);
    return res;
}


