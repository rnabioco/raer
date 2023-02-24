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


regidx_t *regidx_load_payload(char** chroms, int* pos, int* strand,
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

regidx_t *regidx_load_simple(char** chroms, int* start, int* end, int n_sites){

  regidx_t *idx = regidx_init(NULL,NULL,NULL,0,NULL);
  if (!idx) Rf_error("[raer interal] init regidx failed\n");

  char *chr_beg;
  int i, ret;
  for(i = 0; i < n_sites; ++i){
    chr_beg = chroms[i];
    hts_pos_t s = (hts_pos_t) start[i] - 1; // convert to 0 based
    hts_pos_t e = (hts_pos_t) end[i] - 1; // inclusive
    ret = regidx_push(idx, chr_beg, chr_beg + strlen(chr_beg) - 1, s, e, NULL);
    if(ret < 0) Rf_error("[raer internal] index push failed\n");
  }
  return idx;
}

regidx_t *parse_bed_gr(SEXP lst){
  if(Rf_length(lst) != 3){
    Rf_error("'lst' must contain seqnames, start, and end");
  }
  SEXP seqnames = VECTOR_ELT(lst, 0);
  SEXP start = VECTOR_ELT(lst, 1);
  SEXP end = VECTOR_ELT(lst, 2);

  int n = Rf_length(seqnames);
  if (!IS_CHARACTER(seqnames) || n == 0)
    Rf_error("'seqnames' must be character");
  if (!IS_INTEGER(start) || LENGTH(start) != n)
    Rf_error("'start' must be integer of length %d", n);
  if (!IS_INTEGER(end) || LENGTH(end) != n)
    Rf_error("'end' must be integer of length %d", n);

  char **seqnms;
  int i;
  seqnms = (char **) R_alloc(sizeof(const char *), n);
  for (i = 0; i < n; ++i){
    seqnms[i] = (char *) translateChar(STRING_ELT(seqnames, i));
  }

  regidx_t* idx = regidx_load_simple(seqnms, INTEGER(start), INTEGER(end), n);
  return idx;

}

regidx_t *parse_variant_gr(SEXP lst){
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
    Rf_error("'pos' must be integer of length %d", n);
  if (!IS_INTEGER(strand) || LENGTH(pos) != n)
    Rf_error("'strand' must be integer(1 = +, 2 = -) of length %d", n);
  if (!IS_CHARACTER(ref) || LENGTH(pos) != n)
    Rf_error("'ref' must be character of length %d", n);
  if (!IS_CHARACTER(alt) || LENGTH(pos) != n)
    Rf_error("'alt' must be character of length %d", n);
  if (!IS_INTEGER(rowidx) || LENGTH(rowidx) != n)
    Rf_error("'rowidx' must be integer of length %d", n);

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

  regidx_t* idx = regidx_load_payload(seqnms, INTEGER(pos), INTEGER(strand),
                                      r, a, INTEGER(rowidx), n);
  return idx;

}
// make flexible to parse simple 3 column bed  vcf-like structure
regidx_t *regidx_build(SEXP lst, int tbl)
{
  regidx_t* idx;
  if(tbl == 1) {
    idx = parse_variant_gr(lst);
  } else if (tbl == 2) {
    idx = parse_bed_gr(lst);
  } else {
    Rf_error("[raer internal] incorrect tbl specification");
  }
  if (!idx){
    Rf_error("[raer internal] regidx indexing failed");
  }
  return idx;
}
