#include "ullr.h"

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include "bri_index.h"
  #include "bri_get.h"
}

// [[Rcpp::export(rng = false)]]
int build_index(std::string bampath,
                std::string idxpath) {
  const char* cbampath = bampath.c_str();
  const char* cidxpath = idxpath.c_str();
  bam_read_idx_build(cbampath, cidxpath) ;
  return 0;
}

// [[Rcpp::export(rng = false)]]
int fetch_cb_reads(std::string bampath,
                   std::string outpath,
                   std::vector<std::string> vals) {

  // need to refactor input bam_read_idx_get_main
  // to pass more sensible data structure
  std::vector<char*> cstrings{};
  std::string prog("get");
  cstrings.push_back(&prog.front()) ;
  cstrings.push_back(&bampath.front()) ;
  cstrings.push_back(&outpath.front()) ;

  for(auto& string : vals)
    cstrings.push_back(&string.front());

  for(auto i:cstrings){
    Rcpp::Rcout << *i << '\n';
  }

  int ret = 0;
  ret = bam_read_idx_get_main(cstrings.size(), cstrings.data()) ;
  return ret;
}