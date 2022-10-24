#ifndef raer_PDATA_H
#define raer_PDATA_H

#include "Rinternals.h"

#ifdef __cplusplus
extern "C" {
#endif

/* largely templated from Rsamtools approach (scanBam) */


/* arrays of pileup data
   to be stored for each file */
typedef struct  {
  int *pos, *nref, *nvar, *na, *nt, *nc, *ng, *nn, *nx;
  char **seqnames, **strand, **ref, **var;
} _PLP_VECS, *PLP_VECS;


typedef struct {
  int BLOCKSIZE; /* size to grow vectors */
  PLP_VECS pdat;  /* structs to grow dynamically*/
  int icnt, ncnt, nfiles;
  FILE **fps;
  SEXP result;  /* list to return to R, will be populated at end of pileup */
} _PLP_DATA, *PLP_DATA;

SEXP pileup_result_init(int n);

SEXP pileup_template();

PLP_DATA init_PLP_DATA(SEXP result, int n);

int grow_PLP_DATA(PLP_DATA pd, int len);

SEXP get_or_grow_PLP_DATA(PLP_DATA pd, int len);

void finish_PLP_DATA(PLP_DATA pd);

enum {
  SEQNAME_IDX = 0, POS_IDX, STRAND_IDX, REF_IDX, VAR_IDX, NREF_IDX, NVAR_IDX,
  NA_IDX, NT_IDX, NC_IDX, NG_IDX, NN_IDX, NX_IDX
};

/* From Rsamtools
 robust memory re-allocation */

#define _Rs_Realloc(p, n, t)	(t *) _Rs_Realloc_impl(p, n, sizeof(t))

void *_Rs_Realloc_impl(void *p, size_t n, size_t t);

#ifdef __cplusplus
}
#endif

#endif
