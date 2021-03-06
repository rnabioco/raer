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
  for (size_t i = 0; i < cbs.size(); ++i)
    cstrings.push_back(cbs[i].c_str());

  int ret = 0;
  ret = bam_read_idx_get(bampath.c_str(),
                         outpath.c_str(),
                         cstrings.data(),
                         cstrings.size()) ;
  return ret;
}

// [[Rcpp::export(rng = false)]]
DataFrame c_show_index(std::string bamfn, std::string idxfn) {

  bam_read_idx* bri = bam_read_idx_load(bamfn.c_str(), idxfn.c_str());
  if(!bri || bri->record_count < 1){
   stop("issue finding tags for index");
  }
  size_t n = bri->record_count;
  CharacterVector tags(n);
  IntegerVector nt(n);
  for(size_t i = 0; i < n; ++i) {
   tags[i] = bri->records[i].read_name.ptr;
   nt[i] = bri->records[i].n_aln;
  }

  DataFrame res;
  res = DataFrame::create(Named("tag") = tags,
                     Named("n") = nt);

  bam_read_idx_destroy(bri);
  return res;
}

// [[Rcpp::export(rng = false)]]
IntegerMatrix cpp_fill_sparse_matrix(const std::vector<std::vector<int> >& vals,
                                     const std::vector<std::vector<int> >& hits) {

  size_t nv = vals.size();
  size_t nh = hits.size();
  if(nv != nh){
    stop("mismatched vals and hits lists");
  }
  if(nv == 0){
    stop("No entries to populate matrix");
  }
  // determine total # of entries
  size_t n = 0;
  for(size_t i = 0; i < nv; ++i){
    n += hits[i].size();
  }

  IntegerMatrix sm(n, 3);
  size_t pos = 0;
  for(size_t i = 0; i < nv; ++i){
    int vs = vals[i].size();
    int hs = hits[i].size();
    if(vs != hs){
      stop("mismatched vals and hits lists");
    }
    for (size_t j = 0; j < vs; ++j){
      sm(pos, 0) = hits[i][j];
      sm(pos, 1) = i + 1;
      sm(pos, 2) = vals[i][j];
      pos += 1;
    }
  }

  return sm;
}
