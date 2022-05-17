#include <stdlib.h>             /* free */
#include "bedfile.h"
#include "bedidx.h"
// Largely templated on approaches used in Rsamtools
// to load, protect, and reuse fasta/bam/tabix indexes from R

static void _bed_close(void * bed)
{
  bed_destroy(bed);
}

static void _bedfile_close(SEXP ext)
{
  _BED_FILE *ffile = BEDFILE(ext);
  if (ffile->index != NULL)
    _bed_close(ffile->index);
  ffile->index = NULL;
}

static void _bedfile_finalizer(SEXP ext)
{
  if (R_ExternalPtrAddr(ext) == NULL)
    return;
  _bedfile_close(ext);
  _BED_FILE *ffile = BEDFILE(ext);
  R_Free(ffile);
  R_SetExternalPtrAddr(ext, NULL);
}

SEXP bedfile_open(SEXP filename)
{
  if (!IS_CHARACTER(filename) || LENGTH(filename) != 1)
    Rf_error("'file' must be character(1)");

  _BED_FILE* ffile = (_BED_FILE*) R_Calloc(1, _BED_FILE);
  const char *fn = translateChar(STRING_ELT(filename, 0));
  SEXP ext = PROTECT(R_MakeExternalPtr(ffile, R_NilValue, filename));
  R_RegisterCFinalizerEx(ext, _bedfile_finalizer, TRUE);
  UNPROTECT(1);

  ffile->index = bed_read(fn);

  if (ffile->index == NULL) {
    R_Free(ffile);
    Rf_error("'indexBed' indexing failed");
  }

  return ext;
}

SEXP bedfile_close(SEXP ext)
{
  _bedfile_close(ext);
  return ext;
}




// debug function to check loading of bed index

//typedef struct {
//  int n, m;
//  uint64_t *a;
//  int *idx;
//  int filter;
//} bed_reglist_t;
//
//
//#include "htslib/khash.h"
//KHASH_MAP_INIT_STR(reg, bed_reglist_t)
//
//typedef kh_reg_t reghash_t;
//
//SEXP print_bed(SEXP ext){
//  Rprintf("hello");
//  Rprintf("%d", Rf_isNull(ext));
//  if(R_ExternalPtrAddr(ext) == NULL){
//    Rf_error("ptr no addre");
//  }
//  if (BEDFILE(ext) == NULL) {
//    Rf_error("corrupt ext");
//    Rf_PrintValue(ext);
//  }
//
// _BED_FILE *ffile = BEDFILE(ext) ;
//  if (ffile->index == NULL)
//    Rf_error("corrupt index");
// Rprintf("bed_read output %p:", ffile );
//
// reghash_t *h = (reghash_t *)ffile->index;
// bed_reglist_t *p;
// khint_t k;
// int i;
// const char *reg;
// uint32_t beg, end;
//
// if (!h) {
//   Rprintf("Hash table is empty!\n");
//   return ScalarLogical(FALSE);
// }
// for (k = kh_begin(h); k < kh_end(h); k++) {
//   if (kh_exist(h,k)) {
//     reg = kh_key(h,k);
//     Rprintf("Region: '%s'\n", reg);
//     if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
//       Rprintf("Filter: %d\n", p->filter);
//       for (i=0; i<p->n; i++) {
//         beg = (uint32_t)(p->a[i]>>32);
//         end = (uint32_t)(p->a[i]);
//
//         Rprintf("\tinterval[%d]: %d-%d\n",i,beg,end);
//       }
//     } else {
//       Rprintf("Region '%s' has no intervals!\n", reg);
//       return ScalarLogical(FALSE);
//     }
//   }
// }
//
//  return ScalarLogical(TRUE);
//}
