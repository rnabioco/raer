#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include "bri/bri_index.h"
  #include "bri/bri_get.h"
}

// [[Rcpp::export(rng = false)]]
int c_build_index(std::string bampath,
                  std::string idxpath,
                  std::string tag) {
  const char* cbampath = bampath.c_str();
  const char* cidxpath = idxpath.c_str();
  const char* ctag = tag.c_str();
  bam_read_idx_build(cbampath, cidxpath, ctag) ;
  return 0;
}

// [[Rcpp::export(rng = false)]]
int fetch_cb_reads(std::string bampath,
                   std::string outpath,
                   std::vector<std::string> cbs) {

  std::vector<const char*> cstrings;
  for (int i = 0; i < cbs.size(); ++i)
    cstrings.push_back(cbs[i].c_str());

  int ret = 0;
  ret = bam_read_idx_get(bampath.c_str(),
                         outpath.c_str(),
                         cstrings.data(),
                         cstrings.size()) ;
  return ret;
}
