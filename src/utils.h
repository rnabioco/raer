#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <Rinternals.h>
#include "htslib/sam.h"
// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer);

int dist_to_splice(bam1_t* b, int pos, int dist);
int dist_to_indel(bam1_t* b, int pos, int dist);
int trim_pos(bam1_t* b, int pos, int dist);
#endif
