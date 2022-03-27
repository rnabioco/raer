#ifndef raer_pileup_H
#define raer_pileup_H

int run_cpileup(const char** cbampath,
                int n,
                const char* fapath,
                const char* cregion,
                const char* outfn,
                const char* bedfn,
                int* min_reads,
                int max_depth,
                int min_baseQ,
                int min_mapQ,
                int libtype,
                const char* r_flags,
                const char* f_flags,
                int n_align,
                const char* n_align_tag,
                int nmer,
                SEXP ext);

#endif
