#ifndef raer__raer_H
#define raer__raer_H

#include <algorithm>

#include <Rcpp.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"

using namespace Rcpp ;

// class for handling bam file opening and closing
class BamReader {
public:
  samFile* in;
  hts_idx_t* idx;
  BGZF* bz;
  bam_hdr_t* header ;
  BamReader(const std::string& bampath,
            bool check_idx = true,
            int cache_size=10*BGZF_MAX_BLOCK_SIZE) ;

  ~BamReader(){
    hts_idx_destroy(idx);
    sam_close(in);
  }
};

#endif
