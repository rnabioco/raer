#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <htslib/sam.h>
#include <Rinternals.h>

SEXP get_region(SEXP region);

int check_simple_repeat(char** ref, hts_pos_t* ref_len, int pos, int nmer);
int dist_to_splice(bam1_t* b, int pos, int dist);
int dist_to_indel(bam1_t* b, int pos, int dist);
int check_variant_pos(bam1_t* b, int pos, int dist_5p, int dist_3p);
int check_variant_fpos(bam1_t* b, int pos, double fdist_5p, double fdist_3p);
int query_start(bam1_t *b);
int query_end(bam1_t *b);
int read_base_quality(bam1_t* b, float pc, int mq);
int invert_read_orientation(bam1_t* b, int libtype);
int check_splice_overhang(bam1_t* b, int pos, int dist);
double calc_sor(int fwd_ref, int rev_ref, int fwd_alt, int rev_alt);

char *reverse(char *str);
char *get_read(const bam1_t *rec);
char *get_aux_ztag(bam1_t *b, const char tag[2]);
extern unsigned char comp_base[256];

#endif
