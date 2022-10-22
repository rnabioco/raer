#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <Rinternals.h>
#include <htslib/sam.h>
#include <htslib/khash.h>

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

KHASH_SET_INIT_STR(str)
typedef khash_t(str) *strhash_t;
int populate_lookup_from_file(strhash_t lookup, char *fn);

KHASH_SET_INIT_STR(varhash)
typedef khash_t(varhash) *varhash_t;

void clear_varhash_set(varhash_t vhash);
int parse_mismatches(bam1_t* b,const int pos, int n_types, int n_mis);
#endif
