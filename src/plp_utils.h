#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <Rinternals.h>
#include <htslib/sam.h>

// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer);

int check_simple_repeat(char** ref, hts_pos_t* ref_len, int pos, int nmer);
int dist_to_splice(bam1_t* b, int pos, int dist);
int dist_to_indel(bam1_t* b, int pos, int dist);
int trim_pos(bam1_t* b, int pos, int dist_5p, int dist_3p);
int query_start(bam1_t *b);
int query_end(bam1_t *b);
int read_base_quality(bam1_t* b, float pc, int mq);
int invert_read_orientation(bam1_t* b, int libtype);
int check_splice_overhang(bam1_t* b, int pos, int dist);

char *reverse(char *str);
char *get_read(const bam1_t *rec);
char *get_aux_ztag(bam1_t *b, const char tag[2]);
unsigned char comp_base[256];

#endif
