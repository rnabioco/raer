#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "bedfile.h"
#include "utils.h"
/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _raer_c_build_index(SEXP, SEXP);
extern SEXP _raer_fetch_cb_reads(SEXP, SEXP, SEXP);
extern SEXP _raer_get_region(SEXP);
extern SEXP _raer_read_bam(SEXP, SEXP, SEXP, SEXP);
extern SEXP _raer_read_bam_tags(SEXP, SEXP, SEXP, SEXP);
extern SEXP _raer_run_pileup(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bedfile_open(SEXP);
extern SEXP bedfile_close(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_raer_c_build_index",  (DL_FUNC) &_raer_c_build_index,  2},
  {"_raer_fetch_cb_reads", (DL_FUNC) &_raer_fetch_cb_reads, 3},
  {"_raer_get_region",     (DL_FUNC) &_raer_get_region,     1},
  {"_raer_read_bam",       (DL_FUNC) &_raer_read_bam,       4},
  {"_raer_read_bam_tags",  (DL_FUNC) &_raer_read_bam_tags,  4},
  {"_raer_run_pileup",     (DL_FUNC) &_raer_run_pileup,     16},
  {".bedfile_open", (DL_FUNC) &bedfile_open, 1},
  {".bedfile_close", (DL_FUNC) &bedfile_close, 1},
  {".isnull", (DL_FUNC) &isnull, 1},
  {NULL, NULL, 0}
};

void R_init_raer(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
