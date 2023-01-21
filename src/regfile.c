#include <stdlib.h>             /* free */
#include <htslib/regidx.h>
#include "regfile.h"


static inline void free_regidx(void *payload){
  payload_t *pld = *((payload_t**)payload);
  if(pld->alt) free(pld->alt);
  if(pld->ref) free(pld->ref);
  free(pld);
}

void load_payload(payload_t *pld, int strand, char* ref,
                  char* alt, int rowidx){
  pld->strand = strand;
  pld->alt = strdup(alt);
  pld->ref = strdup(ref);
  pld->idx = rowidx;
}

regidx_t *regidx_load(char** chroms, int* pos, int* strand,
                      char** ref, char** alt, int* rowidx,
                      int n_sites){

  regidx_t *idx = regidx_init(NULL,NULL,free_regidx,sizeof(payload_t*),NULL);
  if (!idx) Rf_error("[raer interal] init regidx failed\n");

  char *chr_beg;
  int i, ret;
  payload_t *pld;
  for(i = 0; i < n_sites; ++i){
    chr_beg = chroms[i];
    pld = (payload_t*) R_Calloc(1, payload_t);
    load_payload(pld, strand[i], ref[i], alt[i], rowidx[i]);

    hts_pos_t p = (hts_pos_t) pos[i] - 1; // convert 1 to 0 based
    ret = regidx_push(idx, chr_beg, chr_beg + strlen(chr_beg) - 1, p, p, &pld);
    if(ret < 0) Rf_error("[raer internal] index push failed\n");
  }
  return idx;
}

regidx_t *regidx_build(SEXP lst)
{
  if(Rf_length(lst) != 6){
    Rf_error("'lst' must contain seqnames, pos, strand, ref, alt, and rowidx");
  }
  SEXP seqnames = VECTOR_ELT(lst, 0);
  SEXP pos = VECTOR_ELT(lst, 1);
  SEXP strand = VECTOR_ELT(lst, 2);
  SEXP ref = VECTOR_ELT(lst, 3);
  SEXP alt = VECTOR_ELT(lst, 4);
  SEXP rowidx = VECTOR_ELT(lst, 5);

  int n = Rf_length(seqnames);
  if (!IS_CHARACTER(seqnames) || n == 0)
    Rf_error("'seqnames' must be character");
  if (!IS_INTEGER(pos) || LENGTH(pos) != n)
    Rf_error("'pos' must be integer of length ", n);
  if (!IS_INTEGER(strand) || LENGTH(pos) != n)
    Rf_error("'strand' must be integer(1 = +, 2 = -) of length ", n);
  if (!IS_CHARACTER(ref) || LENGTH(pos) != n)
    Rf_error("'ref' must be character of length ", n);
  if (!IS_CHARACTER(alt) || LENGTH(pos) != n)
    Rf_error("'alt' must be character of length ", n);
  if (!IS_INTEGER(rowidx) || LENGTH(rowidx) != n)
    Rf_error("'rowidx' must be integer of length ", n);

  char **seqnms, **r, **a;
  int i;
  seqnms = (char **) R_alloc(sizeof(const char *), n);
  r = (char **) R_alloc(sizeof(const char *), n);
  a = (char **) R_alloc(sizeof(const char *), n);

  for (i = 0; i < n; ++i){
    seqnms[i] = (char *) translateChar(STRING_ELT(seqnames, i));
    r[i] = (char *) translateChar(STRING_ELT(ref, i));
    a[i] = (char *) translateChar(STRING_ELT(alt, i));
  }

  regidx_t* idx = regidx_load(seqnms, INTEGER(pos), INTEGER(strand),
                 r, a, INTEGER(rowidx), n);

  if (!idx){
    Rf_error("'regidx_build_internal' indexing failed");
  }
  return idx;
}

// Largely templated on approaches used in Rsamtools
// to load, protect, and reuse indexes built in C from R

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

