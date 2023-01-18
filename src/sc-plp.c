#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/hts_log.h>
#include <htslib/regidx.h>
#include <htslib/khash_str2int.h>
#include <Rinternals.h>
#include "regfile.h"
#include "plp.h"
#include "plp_utils.h"



// data structure for storing cell barcodes and UMI information per site
// cbmap:
//   - hashmap with cellbarcode as key,
//   - value is a hashset populated with UMIs as keys

KHASH_SET_INIT_STR(str)
typedef khash_t(str) *strhash_t;

KHASH_MAP_INIT_STR(str2intmap, int)
typedef khash_t(str2intmap) *str2intmap_t;

typedef struct  {
  strhash_t umi; // hashset mapping UMI to consensus calls
  int idx; // cb column index
  int ref; // # of UMIs supporting ref
  int alt; // # of UMIs supporting alt
} cb_t;

KHASH_MAP_INIT_STR(cbumi_map, cb_t)
  typedef khash_t(cbumi_map) *cbumi_map_t;

void clear_umi(strhash_t uhash) {
  khint_t k;
  if(!uhash) return;
  for (k = kh_begin(uhash); k < kh_end(uhash); ++k){
    if (kh_exist(uhash, k)) free((char*)kh_key(uhash, k));
  }
  kh_clear(str, uhash);
}

void clear_cb_umiset(cbumi_map_t cbhash) {
  khint_t k;
  cb_t cdat;
  if(!cbhash) return;
  for (k = kh_begin(cbhash); k < kh_end(cbhash); ++k){
    if (!kh_exist(cbhash, k)) continue;
    cdat = kh_value(cbhash, k);
    clear_umi(cdat.umi);
    cdat.ref = cdat.alt = 0;
  }
}




//From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
static void chkIntFn(void *dummy) {
    R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
int _checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

// struct for event filter params
typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm, min_overhang, min_var_reads;
} efilter;

typedef struct {
  int minq;
  double pct;
} read_qual_t;

typedef struct {
  int min_mq, libtype, min_bq, max_depth;
  read_qual_t read_qual;
  uint32_t keep_flag[2];
  char* qregion;
  regidx_t *reg_idx;
  regitr_t *reg_itr;
  char* mtxfn;
  char* bcfn;
  char* sitesfn;
  int nbam;
  char **bamfns;
  FILE **fps; // [0] = mtxfn, [1] = sitefn, [2] = bcfn;
  efilter ef;
  str2intmap_t cbidx; // hashmap cellbarcode key -> column index value
  cbumi_map_t cbmap; // hashmap cellbarcode key -> cbumi_map_t as value
  int has_umi;
  char* umi_tag;
  char* cb_tag;
  int idx_skip;
  int pe; // 1 if paired end, 0 if not
  int report_all; // 1 report all sites in sparseMatric
  int site_idx; // counter of sites writen, used for row index if report_all is FALSE
} sc_mplp_conf_t;

typedef struct {
  char *ref[2];
  int ref_id[2];
  hts_pos_t ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
  samFile *fp;
  hts_itr_t *iter;
  bam_hdr_t *h;
  mplp_ref_t *ref;
  const sc_mplp_conf_t *conf;
  hts_idx_t *idx;
} mplp_aux_t;

typedef struct {
  int n;
  int *n_plp, *m_plp;
  bam_pileup1_t **plp;
} mplp_pileup_t;

// structs for holding counts
typedef struct {
  int nr, nv, no;
  int ref_b;
} counts;

typedef struct  {
  int pos;
  counts *pc;
  counts *mc;
} pcounts;


/* return codes,
 * 0 = read passes,
 * 1 = read fails filter should be counted as bad read
 * 2 = read fails do no count as bad read
 */
static int check_read_filters(const bam_pileup1_t *p, sc_mplp_conf_t *conf){

  // skip indel and ref skip ;
  if(p->is_del || p->is_refskip) return(2) ;

  // filter based on mapq
  if (p->b->core.qual < conf->min_mq) return(1) ;

  // check base quality
  int bq = p->qpos < p->b->core.l_qseq
  ? bam_get_qual(p->b)[p->qpos]
  : 0;
  if (bq < conf->min_bq) {
    if(bq == 0) { // overlapping read pairs get base qualities set to 0 by mplp
      return(2);
    } else {
      return(1);
    }
  }

  // check if pos is within x dist from 5' end of read, qpos is 0-based
  if(trim_pos(p->b, p->qpos, conf->ef.trim_5p_dist, conf->ef.trim_3p_dist)) return(1);

  // check for splice in alignment nearby
  if(conf->ef.splice_dist && dist_to_splice(p->b, p->qpos, conf->ef.splice_dist) >= 0) return(1);

  // check if site in splice overhang and > min_overhang
  if(conf->ef.min_overhang && check_splice_overhang(p->b, p->qpos, conf->ef.min_overhang) > 0) return(1);

  // check if indel event nearby
  if(conf->ef.indel_dist && dist_to_indel(p->b, p->qpos, conf->ef.indel_dist) >= 0) return(1);

  return 0;
}

// -1 on missing cb or umi
// 0 umi is duplicate
// 1 umi is new
static int count_record(bam1_t *b, sc_mplp_conf_t *conf, payload_t *pld,
                        unsigned char base){
  khiter_t k;
  char *cb, *umi, *umi_val;
  cb_t cbdat;
  int uret;
  cb = get_aux_ztag(b, conf->cb_tag);
  if(cb == NULL) return(-1);
  k = kh_get(cbumi_map, conf->cbmap, cb);
  if(k == kh_end(conf->cbmap)) return(-1);
  cbdat = kh_value(conf->cbmap, k);

  umi = get_aux_ztag(b, conf->umi_tag);
  if(umi == NULL) return(-1);
  umi_val = strdup(umi);
  kh_put(str, cbdat.umi, umi_val, &uret);
  if (uret == 0) {
    free(umi_val);
    return(0);
  }

  if(base == (unsigned char) *pld->ref){
    cbdat.ref += 1;
  } else if (base == (unsigned char) *pld->alt){
    cbdat.alt += 1;
  }
  return(1);
}

// write entry to sparse array-like format
// col 1 = row index (site index)
// col 2 = column index (barcode index)
// col 3 = values for ref
// col 4 = values for alt
//
// read into R as a 4 column matrix, then coerce to list of 2 sparseMatrices.
static int write_counts(sc_mplp_conf_t *conf, payload_t *pld, const char* seqname, int pos){
  const char *cb;
  cb_t cbdat;
  int c_idx, n_rec = 0;
  khint_t k, j;
  for (k = kh_begin(conf->cbmap); k != kh_end(conf->cbmap); ++k) {
    if (!kh_exist(conf->cbmap, k)) continue;
    cb = kh_key(conf->cbmap, k);
    cbdat = kh_val(conf->cbmap, k);
    j = kh_get(str2intmap, conf->cbidx, cb);
    if(j != kh_end(conf->cbidx)) {
      c_idx = kh_value(conf->cbidx, j);
    } else {
      return(-1);
    }
    if(conf->report_all){
      fprintf(conf->fps[0], "%d\t%d\t%d\t%d\n", pld->idx, c_idx, cbdat.ref, cbdat.alt);
    } else {
      fprintf(conf->fps[0], "%d\t%d\t%d\t%d\n", conf->site_idx, c_idx, cbdat.ref, cbdat.alt);
      fprintf(conf->fps[1], "%s_%d_%s_%s_%s\n", seqname, pos, pld->strand, pld->ref, pld->alt);
      conf->site_idx += 1;
    }
    n_rec += 1;
  }
  return(n_rec);
}

// read processing function for pileup
static int screadaln(void *data, bam1_t *b) {
  mplp_aux_t *g = (mplp_aux_t *)data;
  int ret, has_cb, skip = 0;
  uint32_t test_flag;
  char* cb;
  do {
    ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->h, b);
    if (ret < 0) break;

    // exclude unmapped reads
    if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) {
      skip = 1;
      continue;
    }

    // exclude read if it has no cb or cb not in requested cbs
    khiter_t k;
    if(g->conf->cbidx && g->conf->cb_tag){
      cb = get_aux_ztag(b, g->conf->cb_tag);
      if(!cb) {skip = 1; continue;}
      k = kh_get(str2intmap, g->conf->cbidx, cb);
      if(k == kh_end(g->conf->cbidx)){
        skip = 1;
        continue;
      }
    }

    // test required and filter flags
    test_flag = (g->conf->keep_flag[0] & ~b->core.flag) |
      (g->conf->keep_flag[1] & b->core.flag);
    if (~test_flag & 2047u){ skip = 1; continue;}

    // test overlap
    if (g->conf->reg_idx) {
      skip = !regidx_overlap(g->conf->reg_idx,
                      sam_hdr_tid2name(g->h, b->core.tid),
                      b->core.pos,
                      bam_endpos(b),
                      g->conf->reg_itr);
      if (skip) continue;
    }

    skip = 0;
    // check mapping quality
    if (b->core.qual < g->conf->min_mq) {
      skip = 1;
    } else if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)){
      skip = 1;
      // check if read quality is > x in at least y% of read
    } else if(g->conf->read_qual.pct &&
      g->conf->read_qual.minq &&
      !read_base_quality(b, g->conf->read_qual.pct, g->conf->read_qual.minq)) {
      skip = 1;
    }
  } while (skip);
  return ret;
}

static void set_event_filters(efilter* ef, int* event_filters){
  ef->trim_5p_dist = event_filters[0];
  ef->trim_3p_dist = event_filters[1];
  ef->splice_dist = event_filters[2];
  ef->indel_dist = event_filters[3];
  ef->nmer = event_filters[4];
  ef->n_mm_type = event_filters[5];
  ef->n_mm = event_filters[6];
  ef->min_overhang = event_filters[7];
  ef->min_var_reads = event_filters[8];
}

static int run_sc_pileup(sc_mplp_conf_t *conf) {

  hts_set_log_level(HTS_LOG_ERROR);

  mplp_aux_t **data;
  int i, tid, *n_plp, ret = 0;
  hts_pos_t pos, beg0 = 0, end0 = INT32_MAX, ref_len;

  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  bam_mplp_t iter;
  bam_hdr_t *h = NULL;
  char *ref;

  data = calloc(conf->nbam, sizeof(mplp_aux_t*));
  plp = calloc(conf->nbam, sizeof(bam_pileup1_t*));
  n_plp = calloc(conf->nbam, sizeof(int));

  // read the header of each file in the list and initialize data
  for (i = 0; i < conf->nbam; ++i) {
    bam_hdr_t *h_tmp;
    data[i] = calloc(1, sizeof(mplp_aux_t));
    data[i]->fp = sam_open(conf->bamfns[i], "rb");
    if ( !data[i]->fp ) {
      Rf_error("failed to open %s: %s\n", conf->bamfns[i], strerror(errno));
    }

    data[i]->conf = conf;
    data[i]->ref = &mp_ref;

    h_tmp = sam_hdr_read(data[i]->fp);
    if ( !h_tmp ) {
      Rf_error("fail to read the header of %s\n", conf->bamfns[i]);
    }

    if(conf->qregion){
      hts_idx_t *idx = NULL;
      idx = sam_index_load(data[i]->fp, conf->bamfns[i]) ;

      if (idx == NULL) {
        Rf_error("fail to load bamfile index for %s\n", conf->bamfns[i]);
      }

      if ( (data[i]->iter=sam_itr_querys(idx, h_tmp, conf->qregion)) == 0) {
        Rf_error("Fail to parse region '%s' with %s\n", conf->qregion, conf->bamfns[i]);
      }

      if(i == 0){
        beg0 = data[i]->iter->beg;
        end0 = data[i]->iter->end;
      }
      hts_idx_destroy(idx);

    } else {
      data[i]->iter = NULL;
    }

    if(i == 0){
      h = data[i]->h = h_tmp;
    } else {
      bam_hdr_destroy(h_tmp);
      data[i]->h = h;
    }
  }

  iter = bam_mplp_init(conf->nbam, screadaln, (void**)data);
  // enable overlap detection
  if(conf->pe) bam_mplp_init_overlaps(iter) ;
  // set max depth
  bam_mplp_set_maxcnt(iter, conf->max_depth);

  if (!iter) {
    REprintf("issue with iterator");
    ret = -1;
    goto fail;
  }
  // initialize output files
  conf->fps[0] = fopen(R_ExpandFileName(conf->mtxfn), "w");
  conf->fps[1] = fopen(R_ExpandFileName(conf->bcfn), "w");
  conf->fps[2] = fopen(R_ExpandFileName(conf->sitesfn), "w");

  for (i = 0; i < 3; ++i) {
    if (conf->fps[i] == NULL) {
      REprintf("Failed to open file outputfiles\n");
      ret = -1;
      goto fail;
    }
  }

  int last_tid = -1;
  int n_iter = 0;
  while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    if (conf->qregion && (pos < beg0 || pos >= end0)) continue; // not in single region requested
    if (tid < 0) break;

    // check user interrupt
    if (n_iter % 262144 == 0) {
      if(_checkInterrupt()) goto fail;
    }
    // ensure position is in requested intervals
    // if so recover payload with information about strand and ref/var bases
    int ol;
    payload_t *pld;
    if (conf->reg_idx && tid >= 0) {
      ol = regidx_overlap(conf->reg_idx, sam_hdr_tid2name(h, tid),
                          pos, pos + 1, conf->reg_itr);
      if(ol){
        pld = regitr_payload(conf->reg_itr, payload_t*);
      } else {
        continue;
      }
    } else {
      continue;
    }

    // reset cb/umi count structure
    clear_cb_umiset(conf->cbmap);

    // iterate through bam files
    for (i = 0; i < conf->nbam; ++i) {
      int j;
      // iterate through reads that overlap position
      for (j = 0; j < n_plp[i]; ++j) {

        const bam_pileup1_t *p = plp[i] + j;

        // get read base
        int c = p->qpos < p->b->core.l_qseq
        ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
        : 'N';

        // store base counts based on library type and r1/r2 status
        int invert = invert_read_orientation(p->b, conf->libtype);
        if(invert < 0){
          REprintf("[internal] invert read orientation failure %i\n", invert);
          ret = -1;
          goto fail;
        }

        // remove bad reads
        int rret = check_read_filters(p, conf);
        if(rret > 0) continue;

        if(invert) c = (char)comp_base[(unsigned char)c];
        int cret = count_record(p->b, conf, pld, c);
        if(cret < 0) {
          ret = -1;
          goto fail;
        }
      }
    }

    // write records
    int wr = write_counts(conf, pld, sam_hdr_tid2name(h, tid), pos);
    if(wr == -1){
      REprintf("[raer internal] error in writing output.");
      ret = -1;
      goto fail;
    }
    ++n_iter;
  }

  fail:
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);

    for (i = 0; i < conf->nbam; ++i) {
      sam_close(data[i]->fp);
      if(data[i]->iter) hts_itr_destroy(data[i]->iter);
      if(data[i]->idx) hts_idx_destroy(data[i]->idx);
      free(data[i]);
      fclose(conf->fps[i]);
    }

    free(data); free(plp); free(n_plp);
    if(conf->fps) R_Free(conf->fps);
    if(conf->reg_itr) regitr_destroy(conf->reg_itr);

  if(ret < 0) Rf_error("error detected during pileup");
  return ret;
}


static int set_sc_mplp_conf(sc_mplp_conf_t *conf, int n_bams, char** bamfns,
                            int n_outfns, char** outfns, char* qregion, regidx_t* idx,
                            int min_mapQ, int min_baseQ, double* read_bqual_filter,
                            int max_depth, int* b_flags, int* event_filters, int libtype,
                            int n_bcs, char** bcs, char* cbtag, char* umi,
                            int idx_skip, int pe, int report_all){
  int ret = 0;
  conf->bamfns = bamfns;
  conf->nbam = n_bams;

  conf->fps = R_Calloc(n_outfns, FILE*);
  if(conf->fps == NULL) Rf_error("[raer internal] unable to alloc");
  conf->mtxfn = outfns[0];
  conf->sitesfn = outfns[1];
  conf->bcfn = outfns[2];

  if(qregion) conf->qregion = qregion;
  if(idx) {
    conf->reg_idx = idx;
    conf->reg_itr = regitr_init(conf->reg_idx);
  };

  conf->min_mq = (min_mapQ < 0) ? 0 : min_mapQ;
  conf->min_bq = (min_baseQ < 0) ?  0 : min_baseQ;
  conf->max_depth = (!max_depth) ?  10000: max_depth;

  if(b_flags){
    conf->keep_flag[0] = b_flags[0];
    conf->keep_flag[1] = b_flags[1];
  }

  if(read_bqual_filter){
    conf->read_qual.pct = read_bqual_filter[0];
    conf->read_qual.minq = (int)read_bqual_filter[1];
  }

  memset(&conf->ef, 0, sizeof(efilter));
  set_event_filters(&(conf->ef), event_filters);
  int hret = 0;
  khint_t k;
  if(n_bcs > 0){
    conf->cbidx = kh_init(str2intmap);
    for(int i = 0; i < n_bcs; ++i){
      char *cb = strdup(bcs[i]);
      k = kh_put(str2intmap, conf->cbidx, cb, &hret);
      if(hret == -1 || hret == 2){
        Rf_error("[raer internal] unable to populate barcode hashmap");
      } else if (hret == 0) {
        free(cb); // was a duplicate cellbarcode
      } else {
        kh_value(conf->cbidx, k) = i;
      }
    }
  } else {
    conf->cbidx = NULL;
  }

  if(umi){
    conf->has_umi = 1;
    conf->umi_tag = umi;
  } else {
    conf->has_umi = 0;
    conf->umi_tag = NULL;
  }
  conf->cb_tag = cbtag,
  conf->idx_skip = idx_skip;
  conf->pe = pe;
  conf->report_all = report_all;
  conf->site_idx = 0;

  return(ret);
}

static void check_sc_plp_args(SEXP bampaths, SEXP qregion, SEXP region_idx,
                              SEXP barcodes, SEXP cbtag, SEXP event_filters,
                              SEXP min_mapQ, SEXP max_depth, SEXP min_baseQ,
                              SEXP read_bqual_filter, SEXP libtype, SEXP b_flags,
                              SEXP outfns, SEXP umi, SEXP index_skip, SEXP pe,
                              SEXP report_all){

  if(!IS_CHARACTER(bampaths) || (LENGTH(bampaths) < 1)){
    Rf_error("'bampaths' must be character");
  }

  if(!IS_CHARACTER(qregion) || (LENGTH(qregion) > 1)){
    Rf_error("'qregion' must be character of length 0 or 1");
  }

  if(Rf_isNull(region_idx)){
    Rf_error("'region_idx' must be non-null");
  }

  if(!IS_CHARACTER(barcodes) || (LENGTH(barcodes) < 1)){
    Rf_error("'barcodes' must be character");
  }

  if((!IS_CHARACTER(cbtag) || (LENGTH(cbtag) > 1))){
     Rf_error("'cbtag' must be character of length 0 or 1");
  }

  if(!IS_INTEGER(event_filters) || (LENGTH(event_filters) != 9)){
    Rf_error("'event_filters' must be integer of length 9");
  }

  if(!IS_INTEGER(min_mapQ) || (LENGTH(min_mapQ) != 1)){
    Rf_error("'min_mapQ' must be integer(1)");
  }

  if(!IS_INTEGER(max_depth) || (LENGTH(max_depth) != 1)){
    Rf_error("'max_depth' must be integer(1)");
  }

  if(!IS_INTEGER(min_baseQ) || (LENGTH(min_baseQ) != 1)){
    Rf_error("'min_baseQ' must be integer(1)");
  }

  if(!IS_NUMERIC(read_bqual_filter) || (LENGTH(read_bqual_filter) != 2)){
    Rf_error("'read_bqual_filter' must be numeric of length 2");
  }

  if(!IS_INTEGER(libtype) || (LENGTH(libtype) != 1)){
    Rf_error("'lib_type' must be integer(1)");
  }

  if(!IS_INTEGER(b_flags) || (LENGTH(b_flags) != 2)){
    Rf_error("'b_flags' must be integer of length 2");
  }

  if((!IS_CHARACTER(outfns) || (LENGTH(outfns) != 3))){
    Rf_error("'outfns' must be character(3)");
  }

  if((!IS_CHARACTER(umi) || (LENGTH(umi) > 1))){
    Rf_error("'umi' must be character of length 0 or 1");
  }

  if(!IS_LOGICAL(index_skip) || (LENGTH(index_skip) != 1)){
    Rf_error("'index_skip' must be logical(1)");
  }

  if(!IS_LOGICAL(pe) || (LENGTH(pe) != 1)){
    Rf_error("'pe' must be logical(1)");
  }

  if(!IS_LOGICAL(report_all) || (LENGTH(report_all) != 1)){
    Rf_error("'report_all' must be logical(1)");
  }

}




SEXP do_run_scpileup(SEXP bampaths, SEXP query_region, SEXP region_idx,
                   SEXP barcodes, SEXP cbtag, SEXP event_filters, SEXP min_mapQ,
                   SEXP max_depth, SEXP min_baseQ, SEXP read_bqual_filter,
                   SEXP libtype, SEXP b_flags, SEXP outfns, SEXP umi,
                   SEXP index_skip, SEXP pe, SEXP report_all) {

  check_sc_plp_args(bampaths, query_region, region_idx,
                    barcodes, cbtag, event_filters, min_mapQ, max_depth,
                    min_baseQ, read_bqual_filter, libtype,
                    b_flags, outfns, umi, index_skip, pe, report_all);

  int i;
  char ** cbampaths;
  int nbams = Rf_length(bampaths);
  cbampaths = (char**) R_alloc(sizeof( char *), nbams);
  for (i = 0; i < nbams; ++i){
    cbampaths[i] = (char *) translateChar(STRING_ELT(bampaths, i));
  }

  char * cq_region = LENGTH(query_region) == 0 ?
                   NULL : (char *) translateChar(STRING_ELT(query_region, 0));

  char** coutfns;
  int nout = Rf_length(outfns);
  coutfns = (char**) R_alloc(sizeof(char *), nout);
  for (i = 0; i < nout; ++i){
    coutfns[i] = (char *) translateChar(STRING_ELT(outfns, i));
  }

  char** bcs;
  int nbcs = Rf_length(barcodes);
  bcs = (char**) R_alloc(sizeof(char *), nbcs);
  for (i = 0; i < nbcs; ++i){
    bcs[i] = (char *) translateChar(STRING_ELT(barcodes, i));
  }

  char * c_cbtag = (char *) translateChar(STRING_ELT(cbtag, 0));

  char * c_umi = LENGTH(umi) == 0 ?
  NULL : (char *) translateChar(STRING_ELT(umi, 0));

  sc_mplp_conf_t ga;
  memset(&ga, 0, sizeof(sc_mplp_conf_t));

  _REG_IDX *ridx = NULL;
  if(!Rf_isNull(region_idx)){
    ridx = REGIDX(region_idx) ;
    if (ridx->index == NULL){
      Rf_error("Failed to load region index");
    }
    ga.reg_idx = ridx->index;
  }

  int ret, res,idx_skip = 0;

  ret = set_sc_mplp_conf(&ga, nbams, cbampaths, nout, coutfns,
                         cq_region, ridx->index, INTEGER(min_mapQ)[0],
                         INTEGER(min_baseQ)[0], REAL(read_bqual_filter),
                         INTEGER(max_depth)[0], INTEGER(b_flags),
                         INTEGER(event_filters), INTEGER(libtype)[0],
                         nbcs, bcs, c_cbtag, c_umi, LOGICAL(index_skip)[0],
                         LOGICAL(pe)[0], LOGICAL(report_all)[0]);

  // write barcodes file

  res = 1;
  res = run_sc_pileup(&ga);

  if(ga.fps) {
    for (i=0; i < nout; ++i) fclose(ga.fps[i]);
    R_Free(ga.fps);
  }
  if(ga.cbidx) kh_destroy(str2intmap, ga.cbidx);
  if(ga.reg_itr) regitr_destroy(ga.reg_itr);
  if(ga.reg_idx) regidx_destroy(ga.reg_idx);

  return ScalarLogical(res) ;
}


