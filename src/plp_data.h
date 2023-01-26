#ifndef raer_PDATA_H
#define raer_PDATA_H

#include "Rinternals.h"

#ifdef __cplusplus
extern "C" {
#endif

/* largely templated from Rsamtools approach (scanBam) */


/* arrays of pileup data
   to be stored for each sample */
typedef struct  {
  int *pos, *nref, *nvar, *na, *nt, *nc, *ng, *nn, *nx;
  char **seqnames, **strand, **ref, **var;
} _PLP_VECS, *PLP_VECS;

/* arrays of data to be stored across all samples
   to become rowData values in output */
typedef struct  {
  double *rpbz, *vdb;
} _SITE_VECS, *SITE_VECS;

typedef struct {
  int BLOCKSIZE; /* size to grow vectors */
  PLP_VECS pdat;  /* structs to grow dynamically for each sample*/
  SITE_VECS sdat; /* struct to dynamically per site, across all samples */
  int icnt, ncnt, nfiles;
  FILE **fps;
  SEXP result;  /* list to return to R, will be populated at end of pileup */
} _PLP_DATA, *PLP_DATA;

SEXP pileup_result_init(int n);

SEXP pileup_template();
SEXP sitedata_template();

PLP_DATA init_PLP_DATA(SEXP result, int n);

int grow_PLP_DATA(PLP_DATA pd, int len);

SEXP get_or_grow_PLP_DATA(PLP_DATA pd, int len, int lst);

void finish_PLP_DATA(PLP_DATA pd);

enum {
  SEQNAME_IDX = 0, POS_IDX, STRAND_IDX, REF_IDX, VAR_IDX, NREF_IDX, NVAR_IDX,
  NA_IDX, NT_IDX, NC_IDX, NG_IDX, NN_IDX, NX_IDX
};
enum {
  SITE_DATA_LST = 0, PLP_DATA_LST
};
/* From Rsamtools
 robust memory re-allocation */

#define _Rs_Realloc(p, n, t)	(t *) _Rs_Realloc_impl(p, n, sizeof(t))

void *_Rs_Realloc_impl(void *p, size_t n, size_t t);

#ifdef __cplusplus
}
#endif

#endif
