#ifndef raer_pileup_H
#define raer_pileup_H


int run_cpileup(const char** cbampath,
                int n,
                char* fapath,
                char* cregion,
                const char** outfn,
                char* bedfn,
                int min_reads,
                int max_depth,
                int min_baseQ,
                int* min_mapQ,
                int* libtype,
                int* b_flags,
                int* event_filters,
                int* only_keep_variants,
                char* reads_fn,
                char* mismatches,
                double* read_bqual_filter,
                SEXP ext);

#endif
