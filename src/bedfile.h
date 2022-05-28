#ifndef BEDFILE_H
#define BEDFILE_H

#include <Rdefines.h>
#include "samtools/bedidx.h"

typedef struct {
  void *index;
} _BED_FILE;

#define BEDFILE(f) ((_BED_FILE *) R_ExternalPtrAddr(f))

SEXP bedfile_open(SEXP filename);
SEXP bedfile_close(SEXP ext);
SEXP print_bed(SEXP ext);
#endif
