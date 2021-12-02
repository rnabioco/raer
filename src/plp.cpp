#include "kentr.h"
#include "htslib/sam.h"

#include <Rcpp.h>
using namespace Rcpp;

typedef struct mplp_fns {
  const char *fname;
  samFile *fp;
  bam_hdr_t *fp_hdr;
  hts_itr_t *iter;
} mplp_fns;

// read level processing for pileup
static int readaln(void *data, bam1_t *b) {
  mplp_fns *g = (mplp_fns*)data;
  int ret;

  while (1) {
    ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->fp_hdr, b);
    if (ret < 0) break;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

    // check for mapping status
    uint8_t* aux_info = bam_aux_get(b, "NH") ;
    int n_errors = 0 ;
    if(!aux_info && n_errors < 5) {
      warning ("Warning: missing tag info in reads, e.g: %s", bam_get_qname(b)) ;
      n_errors += 1 ;
      continue ;
    }
    int64_t mmap_tag = bam_aux2i(aux_info);
    if(mmap_tag > 1) continue ;

    break;
  }

  return ret;
}

static int test_mpileup(mplp_fns *input) {
  bam_mplp_t iter = NULL;
  const bam_pileup1_t *plp[1] = { NULL };

  // init int array, if multiple files, set to n files
  int n_plp[1] = { 0 };
  int tid, pos, n = 0;
  int max_count = 10000;
  iter = bam_mplp_init(1, readaln, (void **) &input);

  // return type void in rhtslib (htslib v 1.7), current v1.14 htslib returns int
  // enable overlap detection
  bam_mplp_init_overlaps(iter) ;

  // set max depth
  bam_mplp_set_maxcnt(iter, max_count);

  if (!iter) {
    perror("bam_plp_init");
  }
  int min_baseQ = 20;
  int n_records = 0;
  while ((n = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    if (tid < 0) break;
    if (n_records > 10) break ;
    if (n_plp[0] > 1){
      Rcout  << input->fp_hdr->target_name[tid] << ' '  << pos+1 << ' ' << n_plp[0]  << '\n';
      int j;
      for (j = 0; j < n_plp[0]; ++j) {
        const bam_pileup1_t *p = plp[0] + j;

        // check base quality
        int bq = p->qpos < p->b->core.l_qseq
          ? bam_get_qual(p->b)[p->qpos]
        : 0;
        if (bq < min_baseQ) continue ;

        // skip indel and ref skip ;
        if(p->is_del || p->is_refskip) continue ;

        Rcout << bam_get_qname(p->b)  << '\n';
        Rcout << seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)] << '\n' ;
        Rcout << bq << '\n';

        n_records += 1 ;
      }
    }
  }
  if (n < 0) {
    bam_mplp_destroy(iter);
    return 0;
  }

  bam_mplp_destroy(iter);
  return 0;

}

// [[Rcpp::export(rng = false)]]
int run_pileup(std::string bampath,
               std::string fapath,
               std::string region = ".") {
  const char* cbampath = bampath.c_str();
  const char* cregion = region.c_str();
  mplp_fns g = { NULL, NULL, NULL };
  g.fname = cbampath;
  g.fp = sam_open(cbampath, "rb");
  g.fp_hdr = sam_hdr_read(g.fp);
  if(region != "."){
    hts_idx_t *idx = NULL;
    idx = sam_index_load( g.fp, g.fname) ;

    if (idx == NULL) {
      stop("fail to load bamfile index for %s\n", g.fname);
    }

    if ( (g.iter=sam_itr_querys(idx, g.fp_hdr, cregion)) == 0) {
      fprintf(stderr, "[E::%s] fail to parse region '%s' with %s\n", __func__, conf->reg, fn[i]);
      exit(EXIT_FAILURE);
    }
    beg0 = g.>iter->beg;
    end0 = g.iter->end, tid0 = data[i]->iter->tid;
    hts_idx_destroy(idx);

  }
  int status ;
  status = test_mpileup(&g) ;

  bam_hdr_destroy(g.fp_hdr);
  sam_close(g.fp);

  return 0;
}
