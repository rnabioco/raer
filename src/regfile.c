#include <stdlib.h>             /* free */
#include <htslib/regidx.h>
#include "regfile.h"

// Largely templated on approaches used in Rsamtools
// to load, protect, and reuse indexes built in C from R

static inline void free_regidx(void *payload){
  char **dat = (char**)payload;
  free(*dat);
}

SEXP print_regidx(SEXP ext){
  if(R_ExternalPtrAddr(ext) == NULL){
       Rf_error("ptr no addre");
  }
  if (REGIDX(ext) == NULL) {
    Rf_error("corrupt ext");
    Rf_PrintValue(ext);
  }

  _REG_IDX *ridx = REGIDX(ext) ;
  if (ridx->index == NULL)
    Rf_error("corrupt index");

  regidx_t *idx = (regidx_t *)ridx->index;
  regitr_t *itr = regitr_init(idx);
  int r = 1;
  while ( regitr_loop(itr) && r ){
    REprintf("chr=%s  beg=%d  end=%d payload=%d\n",
             itr->seq,
             itr->beg+1,
             itr->end+1,
             regitr_payload(itr, int));
  }
  regitr_destroy(itr);
  return ScalarLogical(TRUE);
}

regidx_t *regidx_load(char** chroms, int* pos, int* rowidx, int n_sites){

  regidx_t *idx = regidx_init(NULL,NULL,NULL,sizeof(int),NULL);
  if (!idx) Rf_error("[raer interal] init regidx failed\n");

  char *chr_beg;
  int i, ret;
  for(i = 0; i < n_sites; ++i){
    chr_beg = chroms[i];
    hts_pos_t p = (hts_pos_t) pos[i] - 1; // convert 1 to 0 based
    ret = regidx_push(idx, chr_beg, chr_beg + strlen(chr_beg) - 1, p, p, (void *)&rowidx[i]);
    if(ret < 0) Rf_error("[raer internal] index push failed\n");
  }
  return idx;
}

static void _regidx_close(SEXP ext)
{
  _REG_IDX *ridx = REGIDX(ext);
  if (ridx->index != NULL)
    regidx_destroy(ridx->index);
  ridx->index = NULL;
}

SEXP regidx_close(SEXP ext)
{
  _regidx_close(ext);
  return ext;
}

static void _regfile_finalizer(SEXP ext)
{
  if (R_ExternalPtrAddr(ext) == NULL)
    return;
  _regidx_close(ext);
  _REG_IDX *ridx = REGIDX(ext);
  R_Free(ridx);
  R_SetExternalPtrAddr(ext, NULL);
}

SEXP regidx_build(SEXP lst)
{
  if(Rf_length(lst) != 3){
    Rf_error("'lst' must contain seqnames, pos, and rowidx");
  }
  SEXP seqnames = VECTOR_ELT(lst, 0);
  SEXP pos = VECTOR_ELT(lst, 1);
  SEXP rowidx = VECTOR_ELT(lst, 2);

  int n = Rf_length(seqnames);
  if (!IS_CHARACTER(seqnames) || n == 0)
    Rf_error("'seqnames' must be character");
  if (!IS_INTEGER(pos) || LENGTH(pos) != n)
    Rf_error("'pos' must be integer of length ", n);
  if (!IS_INTEGER(rowidx) || LENGTH(rowidx) != n)
    Rf_error("'rowidx' must be integer of length ", n);

  _REG_IDX* ridx = (_REG_IDX*) R_Calloc(1, _REG_IDX);
  SEXP ext = PROTECT(R_MakeExternalPtr(ridx, R_NilValue, lst));
  R_RegisterCFinalizerEx(ext, _regfile_finalizer, TRUE);
  UNPROTECT(1);

  char ** seqnms;
  int i;
  seqnms = (char **) R_alloc(sizeof(const char *), n);
  for (i = 0; i < n; ++i){
    seqnms[i] = (char *) translateChar(STRING_ELT(seqnames, i));
  }

  ridx->index = regidx_load(seqnms, INTEGER(pos), INTEGER(rowidx), n);

  if (ridx->index == NULL) {
    R_Free(ridx);
    Rf_error("'regidx_build' indexing failed");
  }

  return ext;
}

//////////

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
