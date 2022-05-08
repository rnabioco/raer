#ifndef raer_cpp_utils_H
#define raer_cpp_utils_H

#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <Rcpp.h>

class TabixReader {
public:
  htsFile* in;
  tbx_t* idx;
  BGZF* bz;
  TabixReader(const std::string& tbxpath,
              bool check_idx = true) ;

  ~TabixReader(){
    tbx_destroy(idx);
    hts_close(in);
  }
};

#endif
