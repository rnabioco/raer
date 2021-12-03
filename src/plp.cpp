#include "rino.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
 #include "cplp.h"
}

// [[Rcpp::export(rng = false)]]
int run_pileup(std::string bampath,
               std::string fapath,
               std::string region,
               std::string outfn,
               std::string bedfn) {
  const char* cbampath = bampath.c_str();
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

  int out;
  out = run_cpileup(cbampath, cfapath, cregion, coutfn, cbedfn);
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
