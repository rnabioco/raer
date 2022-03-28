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
               std::string outfn,
               std::string bedfn,
               std::vector<int> min_reads,
               std::vector<int> event_filters,
               std::vector<int> min_mapQ,
               std::vector<int> bam_flags,
               std::vector<int> libtype,
               int max_depth = 10000,
               int min_baseQ = 20,
               int n_align = 0,
               std::string n_align_tag = "NH",
               SEXP ext = R_NilValue) {

  int n_files;
  n_files = bampaths.size();
  std::vector<const char *> cbampaths;
  cbampaths.reserve(bampaths.size());
  for(int i = 0; i < bampaths.size(); ++i){
    cbampaths.push_back(bampaths[i].c_str());
  }

  char* cfapath = &*fapath.begin();
  char* coutfn = &*outfn.begin();
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

  if(min_reads.size() == 0){
    stop("please supply min_reads parameter");
  }

  int n_params = 0;
  n_params = n_files - min_reads.size();

  // set min reads to first parameter if not supplied for each file
  if(n_params != 0){
    for(int i = 0; i < n_params; i++){
      min_reads.push_back(min_reads[0]);
    }
  }

  if(min_mapQ.size() == 0){
    stop("please supply min_mapQ parameter");
  }

  n_params = 0;
  n_params = n_files - min_mapQ.size();

  // set min_mapQ to first parameter if not supplied for each file
  if(n_params != 0){
    for(int i = 0; i < n_params; i++){
      min_mapQ.push_back(min_mapQ[0]);
    }
  }


  if(event_filters.size() != 4){
    stop("event filters must be a vector of 4 positive integers ");
  }

  int out;
  out = run_cpileup(&cbampaths[0],
                    cbampaths.size(),
                    cfapath,
                    cregion,
                    coutfn,
                    cbedfn,
                    &min_reads[0],
                    max_depth,
                    min_baseQ,
                    &min_mapQ[0],
                    &libtype[0],
                    &*bam_flags.begin(),
                    n_align,
                    &*n_align_tag.begin(),
                    &event_filters[0],
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
  // i'm sure there is an easier way than this:
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
