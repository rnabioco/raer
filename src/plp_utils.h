#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/faidx.h>
#include "regfile.h"
//#include "samtools/samtools-ext.h"
#include <Rinternals.h>

SEXP get_region(SEXP region);

int check_simple_repeat(char** ref, hts_pos_t* ref_len, int pos, int nmer);
int dist_to_splice(bam1_t* b, int pos, int dist);
int dist_to_indel(bam1_t* b, int pos, int dist);
int check_variant_pos(bam1_t* b, int pos, int dist_5p, int dist_3p);
int check_variant_fpos(bam1_t* b, int pos, double fdist_5p, double fdist_3p);
int query_start(bam1_t* b);
int query_end(bam1_t* b);
int read_base_quality(bam1_t* b, float pc, int mq);
int invert_read_orientation(bam1_t* b, int libtype);
int check_splice_overhang(bam1_t* b, int pos, int dist);
double calc_sor(int fwd_ref, int rev_ref, int fwd_alt, int rev_alt);
int get_relative_position(const bam_pileup1_t* p, int nbase_positions);

KHASH_SET_INIT_STR(strset)
typedef khash_t(strset)* strset_t;
void clear_str_set(strset_t s);

KHASH_MAP_INIT_STR(str2intmap, int)
typedef khash_t(str2intmap)* str2intmap_t;
void clear_str2int_hashmap(str2intmap_t vhash);

char* get_aux_ztag(bam1_t* b, const char tag[2]);

extern const char nt5_str[5];

// get base as character
#define nt5_char(i) (nt5_str[i])

// get 0,1,2,3 index from 4bit encoded base
#define nt16_idx(i) (seq_nt16_int[i])


extern unsigned char comp_base[256];

typedef struct {
  double f5p;
  double f3p;
  int i5p;
  int i3p;
} trim_t;

typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm, min_overhang, min_var_reads;
} efilter;

typedef struct {
  int minq;
  double pct;
} read_qual_t;

typedef struct {
  int min_global_mq, flag, min_bq, min_depth, max_depth, output_reads;
  int report_multiallelics, multi_itr, in_memory;
  int nmer, splice_dist, indel_dist;
  int n_mm_type, n_mm, min_overhang, min_var_reads;
  int nbam, nfps;
  double min_af;
  int umi;
  char* umi_tag;
  int* min_mqs; // across all bam files
  int* libtype; // across all bam files
  int* only_keep_variants; // across all bam files
  trim_t trim;
  read_qual_t read_qual;
  uint32_t keep_flag[2];
  char* reg, *fai_fname, *output_fname;
  faidx_t* fai;
  regidx_t* reg_idx;
  regitr_t* reg_itr;
} mplp_conf_t;

//From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
void chkIntFn(void* dummy) ;
// this will call the above in a top-level context so it won't longjmp-out of context
int checkInterrupt(void);


#endif
