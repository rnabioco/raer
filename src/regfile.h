#ifndef RAER_REGFILE_H
#define RAER_REGFILE_H

#include <Rdefines.h>
#include <bedidx.h> /* packaged with Rhtslib */
#include <htslib/regidx.h>

typedef struct {
  char *strand, *ref, *alt;
  int idx;
} payload_t;

typedef struct {
  regidx_t *index;
} _REG_IDX;

typedef struct {
  void *index;
} _BED_FILE;


#define REGIDX(f) ((_REG_IDX *) R_ExternalPtrAddr(f))
#define BEDFILE(f) ((_BED_FILE *) R_ExternalPtrAddr(f))

SEXP regidx_build(SEXP lst);
SEXP regidx_close(SEXP ext);
SEXP print_regidx(SEXP ext);

SEXP bedfile_open(SEXP);
SEXP bedfile_close(SEXP);
SEXP print_bed(SEXP);


#endif
