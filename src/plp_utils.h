#ifndef raer_UTILS_H
#define raer_UTILS_H

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/faidx.h>
#include "regfile.h"
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

extern unsigned char comp_base[256];

typedef struct  {
    int nmer, splice_dist, indel_dist, trim_i5p, trim_i3p;
    int n_mm_type, n_mm, min_overhang, min_var_reads;
    double trim_f5p, trim_f3p;
} efilter_t;

typedef struct {
  int minq;
  double pct;
} read_qual_t;

typedef struct {
  int min_global_mq;        // min_mapQ across all bams
  int min_bq;               // min base quality
  int min_depth;            // min read coverage 
  int max_depth;            // max # of reads to consider at a site 
  int report_multiallelics; // if 1, report multiple variants per site, otherwise discard site
  int remove_overlaps;      // if 1, enable overlap detection
  int nbam;                 // # of input bam files
  double min_af;            // min allelic frequency required to report variant
  int umi;                  // if 1, check umi tag
  char* umi_tag;            // bam tag containing UMI 
  int* min_mqs;             // min mapQ values for all bam files
  int* libtype;             // library type for all bam files
  int* only_keep_variants;  // if 1, only report variants, setting applies across all bam files
  efilter_t ef;             // various additional site filters 
  read_qual_t read_qual;    // read level base quality filter
  uint32_t keep_flag[2];    // bam flag filter
  char* reg;                // single region to query 
  char* fai_fname;          // indexed fasta filename
  faidx_t* fai;             // fasta index
  regidx_t* reg_idx;        // multiple region index 
  regitr_t* reg_itr;        // multiple region iterator 
} mplp_conf_t;

//From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
void chkIntFn(void* dummy) ;
// this will call the above in a top-level context so it won't longjmp-out of context
int checkInterrupt(void);


#endif
