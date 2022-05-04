#include "raer.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
 #include "cplp.h"
}

//[[Rcpp::export(rng = false)]]
int run_pileup(std::vector<std::string> bampaths,
               std::string fapath,
               std::string region,
               std::string bedfn,
               int min_reads,
               std::vector<int> event_filters,
               std::vector<int> min_mapQ,
               std::vector<int> bam_flags,
               std::vector<int> libtype,
               std::vector<std::string> outfns,
               int max_depth = 10000,
               int min_baseQ = 20,
               int only_keep_variants = 0,
               std::string reads = ".",
               SEXP ext = R_NilValue) {

  int n_files;
  n_files = bampaths.size();
  std::vector<const char *> cbampaths;
  cbampaths.reserve(bampaths.size());
  for(int i = 0; i < bampaths.size(); ++i){
    cbampaths.push_back(bampaths[i].c_str());
  }

  char* cfapath = &*fapath.begin();

  std::vector<const char *> coutfns;
  coutfns.reserve(outfns.size());
  for(int i = 0; i < outfns.size(); ++i){
    coutfns.push_back(outfns[i].c_str());
  }

  char* cregion;
  if(region == "."){
    cregion = NULL;
  } else {
    cregion = &*region.begin();
  }

  char* cbedfn;
  if(bedfn == "."){
    cbedfn = NULL;
  } else {
    cbedfn = &*bedfn.begin();
  }

  if(min_reads < 0){
    stop("min_reads must be positive");
  }

  if(min_mapQ.size() == 0){
    stop("please supply min_mapQ parameter");
  }

  int n_params = 0;
  n_params = n_files - min_mapQ.size();

  // set min_mapQ to first parameter if not supplied for each file
  if(n_params != 0){
    for(int i = 0; i < n_params; i++){
      min_mapQ.push_back(min_mapQ[0]);
    }
  }


  if(event_filters.size() != 7){
    stop("event filters must be a vector of 7 positive integers ");
  }

  char* creadsoutfn;
  if(reads == "."){
    creadsoutfn = NULL;
  } else {
    creadsoutfn = &*reads.begin();
  }

  int out;
  out = run_cpileup(&cbampaths[0],
                    cbampaths.size(),
                    cfapath,
                    cregion,
                    &coutfns[0],
                    cbedfn,
                    min_reads,
                    max_depth,
                    min_baseQ,
                    &min_mapQ[0],
                    &libtype[0],
                    &*bam_flags.begin(),
                    &event_filters[0],
                    only_keep_variants,
                    creadsoutfn,
                    ext);

  return out;
}

// [[Rcpp::export(rng = false)]]
List get_region(std::string region){
  const char* cregion;
  cregion = region.c_str() ;
  int beg, end;
  const char *chr_pos ;
  chr_pos = hts_parse_reg(cregion, &beg, &end) ;
  if(!chr_pos){
    stop("could not parse region:%s", region);
  }
  char *chr_name = (char*)malloc(chr_pos - cregion + 1);
  memcpy(chr_name, cregion, chr_pos - cregion);
  chr_name[chr_pos - cregion] = '\0';
  String chr_name_r(chr_name);

  List res;
  // return 0 based start
  res = List::create(Named("chrom") = chr_name_r,
                     Named("start") = beg,
                     Named("end") = end);
  return res;
}
