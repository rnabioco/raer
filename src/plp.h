#ifndef raer_pileup_H
#define raer_pileup_H

#include <Rinternals.h>
#include <htslib/sam.h>

typedef struct {
  int n;
  int *n_plp, *m_plp;
  bam_pileup1_t **plp;
} mplp_pileup_t;

typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm, min_overhang, min_var_reads;
} efilter;

typedef struct {
  int minq;
  double pct;
} read_qual_t;

SEXP pileup(char** cbampath,
            int n,
            char* fapath,
            char* cregion,
            int in_mem,
            int multi_region_itr,
            const char** outfn,
            const char* bedfn,
            int min_reads,
            int max_depth,
            int min_baseQ,
            int* min_mapQ,
            int* libtype,
            int* b_flags,
            int* event_filters,
            int* only_keep_variants,
            const char* reads_fn,
            char* mismatches,
            double* read_bqual_filter,
            SEXP ext);

#endif
