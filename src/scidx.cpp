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

// [[Rcpp::export(rng = false)]]
IntegerMatrix cpp_fill_sparse_matrix(const std::vector<std::vector<int> >& vals,
                                     const std::vector<std::vector<int> >& hits) {

  int nv = vals.size();
  int nh = hits.size();
  if(nv != nh){
    stop("mismatched vals and hits lists");
  }
  if(nv == 0){
    stop("No entries to populate matrix");
  }
  // determine total # of entries
  int n = 0;
  for(int i = 0; i < nv; ++i){
    n += hits[i].size();
  }

  IntegerMatrix sm(n, 3);
  int pos = 0;
  for(int i = 0; i < nv; ++i){
    int vs = vals[i].size();
    int hs = hits[i].size();
    if(vs != hs){
      stop("mismatched vals and hits lists");
    }
    for (int j = 0; j < vs; ++j){
      sm(pos, 0) = hits[i][j];
      sm(pos, 1) = i + 1;
      sm(pos, 2) = vals[i][j];
      pos += 1;
    }
  }

  return sm;
}
