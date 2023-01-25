#ifndef RAER_REGFILE_H
#define RAER_REGFILE_H

#include <Rdefines.h>
#include <bedidx.h> /* packaged with Rhtslib */
#include <htslib/regidx.h>

typedef struct {
  char *ref, *alt;
  int idx, strand;
} payload_t;

typedef struct {
  void *index;
} _BED_FILE;

#define BEDFILE(f) ((_BED_FILE *) R_ExternalPtrAddr(f))

regidx_t *regidx_build(SEXP lst, int tbl);

SEXP bedfile_open(SEXP);
SEXP bedfile_close(SEXP);
SEXP print_bed(SEXP);


#endif
