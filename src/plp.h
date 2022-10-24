#ifndef raer_pileup_H
#define raer_pileup_H

#include <Rinternals.h>
SEXP run_cpileup(char** cbampath,
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
