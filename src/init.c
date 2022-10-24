#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "bedfile.h"
#include "plp_utils.h"

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _raer_c_build_index(SEXP, SEXP, SEXP);
extern SEXP _raer_cpp_fill_sparse_matrix(SEXP, SEXP);
extern SEXP _raer_cread_tabix(SEXP, SEXP);
extern SEXP _raer_fetch_cb_reads(SEXP, SEXP, SEXP);
extern SEXP _raer_get_region(SEXP);
extern SEXP _raer_list_tabix_chroms(SEXP);
extern SEXP _raer_c_show_index(SEXP, SEXP);
extern SEXP bedfile_open(SEXP);
extern SEXP bedfile_close(SEXP);
extern SEXP do_run_pileup(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_raer_c_build_index",  (DL_FUNC) &_raer_c_build_index,  3},
  {"_raer_cpp_fill_sparse_matrix", (DL_FUNC) &_raer_cpp_fill_sparse_matrix, 2},
  {"_raer_cread_tabix",   (DL_FUNC) &_raer_cread_tabix,    2},
  {"_raer_fetch_cb_reads", (DL_FUNC) &_raer_fetch_cb_reads, 3},
  {"_raer_get_region",     (DL_FUNC) &_raer_get_region,     1},
  {"_raer_list_tabix_chroms", (DL_FUNC) &_raer_list_tabix_chroms,     1},
  {"_raer_c_show_index", (DL_FUNC) &_raer_c_show_index,     2},
  {".bedfile_open", (DL_FUNC) &bedfile_open, 1},
  {".bedfile_close", (DL_FUNC) &bedfile_close, 1},
  {".isnull", (DL_FUNC) &isnull, 1},
  {".do_run_pileup",(DL_FUNC) &do_run_pileup, 20},
  {NULL, NULL, 0}
};

void R_init_raer(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
