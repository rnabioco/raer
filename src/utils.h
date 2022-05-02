#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <Rinternals.h>
#include "htslib/sam.h"
// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer);

int check_simple_repeat(char** ref, int* ref_len, int pos, int nmer);
int dist_to_splice(bam1_t* b, int pos, int dist);
int dist_to_indel(bam1_t* b, int pos, int dist);
int trim_pos(bam1_t* b, int pos, int dist_5p, int dist_3p);
int query_start(bam1_t *b);
int query_end(bam1_t *b);

// From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
inline static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
inline int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}
#endif
