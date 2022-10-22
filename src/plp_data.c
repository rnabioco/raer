#include <Rinternals.h>
#include "Rdefines.h"
#include "plp_data.h"

/* All of the following code mostly is templated on the Rsamtools approach
   to store and grow a c-level datastructure with bam/pileup data.
   The main deviation is handling data from multiple bam files, and the fields stores.
   The datastructure grows per position across multiple files, in contrast to
   growing per region. On finish, the c-level datastucture will be converted to a SEXP,
   avoiding realloc'ing a SEXP during the pileup.
 */

static const int BAM_INIT_SIZE = 1048576;

PLP_DATA init_PLP_DATA(SEXP result, int n) {
  int i;
  PLP_DATA plpd = R_Calloc(1, _PLP_DATA);
  plpd->pdat = R_Calloc(n, _PLP_VECS);
  for(i = 0; i < n; ++i){
    plpd->pdat[i].na = R_Calloc(1, int);
    plpd->pdat[i].nt = R_Calloc(1, int);
    plpd->pdat[i].ng = R_Calloc(1, int);
    plpd->pdat[i].nc = R_Calloc(1, int);
    plpd->pdat[i].nn = R_Calloc(1, int);
    plpd->pdat[i].nx = R_Calloc(1, int);
    plpd->pdat[i].nref = R_Calloc(1, int);
    plpd->pdat[i].nvar = R_Calloc(1, int);
    plpd->pdat[i].pos = R_Calloc(1, int);
    plpd->pdat[i].ref = R_Calloc(1, char*);
    plpd->pdat[i].var = R_Calloc(1, char*);
    plpd->pdat[i].seqnames = R_Calloc(1, char*);
    plpd->pdat[i].strand = R_Calloc(1, char*);
  }
  plpd->BLOCKSIZE = BAM_INIT_SIZE;
  plpd->result = result;
  plpd->nfiles = n;
  return plpd;
}

/* Realloc c-level datastructure holding pileup counts */
int grow_PLP_DATA(PLP_DATA pd, int len)
{
  int i,j;
  SEXP r, s;

  for(i = 0; i < pd->nfiles; ++i){
    r = VECTOR_ELT(pd->result, i);
    for (j = 0; j < LENGTH(r); ++j) {
      if (R_NilValue == (s = VECTOR_ELT(r, j)))
        continue;
      switch (j) {
      case SEQNAME_IDX:
        pd->pdat[i].seqnames = _Rs_Realloc(pd->pdat[i].seqnames, len, char *);
        break;
      case POS_IDX:
        pd->pdat[i].pos = _Rs_Realloc(pd->pdat[i].pos, len, int);
        break;
      case STRAND_IDX:
        pd->pdat[i].strand = _Rs_Realloc(pd->pdat[i].strand, len, char *);
        break;
      case REF_IDX:
        pd->pdat[i].ref = _Rs_Realloc(pd->pdat[i].ref, len, char *);
        break;
      case VAR_IDX:
        pd->pdat[i].var = _Rs_Realloc(pd->pdat[i].var, len, char *);
        break;
      case NREF_IDX:
        pd->pdat[i].nref = _Rs_Realloc(pd->pdat[i].nref, len, int);
        break;
      case NVAR_IDX:
        pd->pdat[i].nvar = _Rs_Realloc(pd->pdat[i].nvar, len, int);
        break;
      case NA_IDX:
        pd->pdat[i].na = _Rs_Realloc(pd->pdat[i].na, len, int);
        break;
      case NT_IDX:
        pd->pdat[i].nt = _Rs_Realloc(pd->pdat[i].nt, len, int);
        break;
      case NC_IDX:
        pd->pdat[i].nc = _Rs_Realloc(pd->pdat[i].nc, len, int);
        break;
      case NG_IDX:
        pd->pdat[i].ng = _Rs_Realloc(pd->pdat[i].ng, len, int);
        break;
      case NN_IDX:
        pd->pdat[i].nn = _Rs_Realloc(pd->pdat[i].nn, len, int);
        break;
      case NX_IDX:
        pd->pdat[i].nx = _Rs_Realloc(pd->pdat[i].nx, len, int);
        break;
      default:
        Rf_error("[raer internal] unhandled grow_PLP_DATA");
      break;
      }
    }
  }

  return len;
}

/* Run at end of pileup
 * Convert the c-data structure
 * into an R list of lists (to become list of GRanges)
 * and free associated memory
 */
void finish_PLP_DATA(PLP_DATA pd) {

  int f_idx, col_idx, row_idx;
  SEXP r, s;
  for(f_idx = 0; f_idx < pd->nfiles; ++f_idx){
    r = VECTOR_ELT(pd->result, f_idx);

    for (col_idx = 0; col_idx < LENGTH(r); ++col_idx) {
      if (R_NilValue == (s = VECTOR_ELT(r, col_idx)))
        continue;

      switch (col_idx) {
      case SEQNAME_IDX:
        //reset length of s to icnt
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        // convert c char* to R strings
        for (row_idx = 0; row_idx < pd->icnt; ++row_idx) {
          SET_STRING_ELT(s, row_idx, mkChar(pd->pdat[f_idx].seqnames[row_idx]));
          R_Free(pd->pdat[f_idx].seqnames[row_idx]);
        }
        R_Free(pd->pdat[f_idx].seqnames);
        break;

      case POS_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        // copy pos array into s
        memcpy(INTEGER(s), pd->pdat[f_idx].pos, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].pos);
        break;

      case STRAND_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        for (row_idx = 0; row_idx < pd->icnt; ++row_idx) {
          SET_STRING_ELT(s, row_idx, mkChar(pd->pdat[f_idx].strand[row_idx]));
          R_Free(pd->pdat[f_idx].strand[row_idx]);
        }
        R_Free(pd->pdat[f_idx].strand);
        break;

      case REF_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        for (row_idx = 0; row_idx < pd->icnt; ++row_idx) {
          SET_STRING_ELT(s, row_idx, mkChar(pd->pdat[f_idx].ref[row_idx]));
          R_Free(pd->pdat[f_idx].ref[row_idx]);
        }
        R_Free(pd->pdat[f_idx].ref);
        break;

      case VAR_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        for (row_idx = 0; row_idx < pd->icnt; ++row_idx) {
          SET_STRING_ELT(s, row_idx, mkChar(pd->pdat[f_idx].var[row_idx]));
          R_Free(pd->pdat[f_idx].var[row_idx]);
        }
        R_Free(pd->pdat[f_idx].var);
        break;

      case NREF_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nref, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nref);
        break;

      case NVAR_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nvar, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nvar);
        break;

      case NA_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].na, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].na);
        break;

      case NT_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nt, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nt);
        break;

      case NC_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nc, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nc);
        break;

      case NG_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].ng, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].ng);
        break;

      case NN_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nn, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nn);
        break;

      case NX_IDX:
        s = Rf_lengthgets(s, pd->icnt);
        SET_VECTOR_ELT(r, col_idx, s);
        memcpy(INTEGER(s), pd->pdat[f_idx].nx, pd->icnt * sizeof(int));
        R_Free(pd->pdat[f_idx].nx);
        break;


      default:
        Rf_error("[raer internal] unhandled finish_PLP_DATA");
      break;
      }
    }
  }
  pd->icnt = pd->ncnt = 0;
}

/* determine if plp_data needs to be realloc
   if not, then return first list element
   list element will be checked for nullness,
   which provides mechanism to select what values should be stored
   For now, storing all specified,
   but could in the future be set in pileup_result_init  */
SEXP get_or_grow_PLP_DATA(PLP_DATA pd, int len)
{
  if (len < 0) {
    if (pd->icnt < pd->ncnt)
      return VECTOR_ELT(pd->result, 0);
    len = pd->ncnt + pd->BLOCKSIZE;
  }

  pd->ncnt = grow_PLP_DATA(pd, len);
  return VECTOR_ELT(pd->result, 0);
}

/* names matching enum defined in plp_data.h
 * will be used as output column names
 */
static const char *TMPL_ELT_NMS[] = {
  "seqname", "pos", "strand", "Ref", "Var", "nRef", "nVar", "nA",
  "nT", "nC", "nG", "nN", "nX"
  /* "vtype", "value" */
};

static const int N_TMPL_ELTS = sizeof(TMPL_ELT_NMS) / sizeof(const char *);

/* Init a list of lists,
   storing a list of pileup vectors from each bamfile.
   Each element of the outer list is data from each bamfile
   The inner list elements contains vectors of pileup data
   On return to R will be list of lists coerced into list of GRanges
 */
SEXP pileup_result_init(int n){

  SEXP result = PROTECT(NEW_LIST(n));
  int i;
  for (i = 0; i < n; ++i) {
    SEXP tmpl = PROTECT(pileup_template());
    SET_VECTOR_ELT(result, i, tmpl);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return result;
}


/* inner list template */
SEXP pileup_template() {

 SEXP tmpl = PROTECT(NEW_LIST(N_TMPL_ELTS));
 SET_VECTOR_ELT(tmpl, SEQNAME_IDX, NEW_CHARACTER(0));
 SET_VECTOR_ELT(tmpl, POS_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, STRAND_IDX, NEW_CHARACTER(0));
 SET_VECTOR_ELT(tmpl, REF_IDX, NEW_CHARACTER(0));
 SET_VECTOR_ELT(tmpl, VAR_IDX, NEW_CHARACTER(0));
 SET_VECTOR_ELT(tmpl, NREF_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NVAR_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NA_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NT_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NC_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NG_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NN_IDX, NEW_INTEGER(0));
 SET_VECTOR_ELT(tmpl, NX_IDX, NEW_INTEGER(0));

 SEXP names = PROTECT(NEW_CHARACTER(N_TMPL_ELTS));
 for (int i = 0; i < N_TMPL_ELTS; ++i)
   SET_STRING_ELT(names, i, mkChar(TMPL_ELT_NMS[i]));
 SET_ATTR(tmpl, R_NamesSymbol, names);
 UNPROTECT(2);
 return tmpl;
}

/* from Rsamtools */
void *_Rs_Realloc_impl(void *p, size_t n, size_t t)
{
  /* Realloc(p, 0, *) fails inappropriately */
  if (n == 0) {
    R_Free(p);
    p = NULL;
  } else {
    p = R_chk_realloc((void *) p, (size_t) (n * t));
  }
  return p;
}

