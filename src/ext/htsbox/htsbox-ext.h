#ifndef RAER_HTSBOX_EXT_H
#define RAER_HTSBOX_EXT_H

#include <htslib/sam.h>
int parse_mismatches(bam1_t* b, const int pos, int n_types, int n_mis);

#endif
