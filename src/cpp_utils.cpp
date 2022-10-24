#include "cpp_utils.h"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include <Rcpp.h>
using namespace Rcpp;


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
  free(chr_name);
  return res;
}


std::vector<std::string> tsv_values(std::string tsv) {

  std::vector<std::string> values ;
  std::stringstream ss(tsv) ;

  while (ss.good()) {
    std::string substr ;
    getline(ss, substr, '\t') ;

    if (substr.empty()) break ;

    values.push_back(substr) ;
  }

  return values ;
}

TabixReader::TabixReader(const std::string& tbxpath,
                         bool check_idx){
  const char* cpath = tbxpath.c_str();

  in = hts_open(cpath, "rz");
  if (in == NULL) {
    stop("Failed to open tabix file, check filepath " + tbxpath);
  }

  bz = in->fp.bgzf ; // bgzf file pointer
  idx = tbx_index_load(cpath); // load tbx index

  if (idx == 0 && check_idx) {
    stop("tabix index file is not available for " + tbxpath);
  }
}

// [[Rcpp::export(rng = false)]]
List cread_tabix(std::string tbxpath,
                std::string region = "."){
  // open tabix
  TabixReader tfile(tbxpath, true) ;

  hts_itr_t *iter = NULL;
  iter = tbx_itr_querys(tfile.idx, region.c_str());
  // initialize empty data record
  kstring_t str = {0,0,0} ;
  int n_fields ;

  std::vector<std::vector<std::string>> output;

  // determine # of output columns
  if (tbx_itr_next(tfile.in, tfile.idx, iter, &str) >= 0){
    std::vector<std::string> fields = tsv_values(str.s) ;
    n_fields = fields.size() ;
    output.resize(n_fields) ;

    for(int i = 0; i < n_fields; i++){
      output[i].push_back(fields[i]) ;
    }

  } else {
    // return empty df
    std::vector<std::string> c;
    std::vector<double> s;
    std::vector<double> e;

    return Rcpp::DataFrame::create(_("chrom") = c ,
                                   _("start") = s,
                                   _("end") = e) ;

  }

  while(tbx_itr_next(tfile.in, tfile.idx, iter, &str) >= 0) {
    std::vector<std::string> fields = tsv_values(str.s) ;
    if(fields.size() != n_fields){
      Rcpp::Rcout << str.s << std::endl ;
      Rcpp::Rcout << n_fields << " " << fields.size() << " "  << std::endl ;
      for(auto i:fields) {
        Rcpp::Rcout << i << std::endl ;
      }
      Rcpp::stop("mismatch field counts") ;
    }

    for(int i = 0; i < n_fields; i++){
      output[i].push_back(fields[i]) ;
    }

  }

  // figure out output names for chrom start end and pos
  int n_no_names = 0 ;
  std::vector<std::string> names(n_fields);
  std::string start_col = "start" ;
  bool end_col = true ;
  if(tfile.idx->conf.bc  == tfile.idx->conf.ec){
    start_col = "pos" ;
    end_col = false ;
  }

  for(int i = 0; i < n_fields; i++){
    if(i == (tfile.idx->conf.sc - 1)){
      names[i] = "chrom" ;
    } else if (i == (tfile.idx->conf.bc - 1)){
      names[i] = start_col;
    } else if (end_col && i == (tfile.idx->conf.ec - 1)) {
      names[i] = "end" ;
    } else {
      n_no_names += 1 ;
      names[i] = "X" + std::to_string(n_no_names) ;
    }
  }

  List res_lst(n_fields);
  for(int i = 0; i < n_fields; i++){
    res_lst[i] = output[i] ;
  }

  DataFrame res = DataFrame::create(res_lst,
                                    _("stringsAsFactors") = false);
  res.attr("names") = names ;

  return res;
}


// [[Rcpp::export(rng = false)]]
CharacterVector list_tabix_chroms(std::string tbxpath){

  TabixReader tfile(tbxpath, true) ;
  int n;
  const char **c_chroms = NULL;
  c_chroms = tbx_seqnames(tfile.idx, &n);
  if (!c_chroms) Rcpp::stop("Failed to get sequence names list");

  CharacterVector chroms(n);

  for (int i = 0; i < n; i++) {
    chroms[i] = c_chroms[i];
  }

  free(c_chroms);
  return chroms;
}
