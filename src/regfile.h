#ifndef RAER_REGFILE_H
#define RAER_REGFILE_H

#include <Rdefines.h>
#include <htslib/regidx.h>

typedef struct {
  char *ref, *alt;
  int idx, strand;
} payload_t;

regidx_t *regidx_build(SEXP lst, int tbl);

#endif
