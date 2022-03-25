#include "raer.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
 #include "cplp.h"
}

// [[Rcpp::export(rng = false)]]
int run_pileup(std::vector<std::string> bampaths,
               std::string fapath,
               std::string region,
               std::string outfn,
               std::string bedfn,
               std::vector<int> min_reads,
               int max_depth = 10000,
               int min_baseQ = 20,
               std::string libtype = "fr-first-strand",
               SEXP ext = R_NilValue) {

  int n_files;
  n_files = bampaths.size();
  std::vector<const char *> cbampaths;
  cbampaths.reserve(bampaths.size());
  for(int i = 0; i < bampaths.size(); ++i){
    cbampaths.push_back(bampaths[i].c_str());
  }

  const char* cfapath = fapath.c_str();
  const char* coutfn = outfn.c_str();
  const char* cregion;
  if(region == "."){
    cregion = NULL;
  } else {
    cregion = region.c_str();
  }

  const char* cbedfn;
  if(bedfn == "."){
    cbedfn = NULL;
  } else {
    cbedfn = bedfn.c_str();
  }

  // encode libtype as 0 = unstranded, 1 = fr-first-strand, 2 = fr-second-strand
  int lib_spec = 0;
  if(libtype == "fr-first-strand"){
    lib_spec = 1;
  } else if (libtype == "fr-second-strand") {
    lib_spec = 2;
  } else if (libtype == "unstranded") {
    lib_spec = 0;
  } else {
    stop("unrecognized library type: fr-first-strand, fr-second-strand, or unstranded supported");
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

  int out;
  out = run_cpileup(&cbampaths[0], cbampaths.size(),
                    cfapath, cregion, coutfn, cbedfn,
                    &min_reads[0], max_depth, min_baseQ, lib_spec, ext);
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
