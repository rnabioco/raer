#ifndef raer_pileup_H
#define raer_pileup_H

int run_cpileup(const char* cbampath,
                const char* fapath,
                const char* cregion,
                const char* outfn,
                const char* bedfn,
                int min_reads,
                int max_depth,
                int min_baseQ,
                int libtype);

#endif
