#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/hts_log.h>
#include <htslib/regidx.h>
#include <htslib/khash_str2int.h>
#include <Rinternals.h>
#include "regfile.h"
#include "plp.h"

KHASH_MAP_INIT_STR(cb2int, int)

  //From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
  static void chkIntFn(void *dummy) {
    R_CheckUserInterrupt();
  }

// this will call the above in a top-level context so it won't longjmp-out of your context
int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

// struct for event filter params
typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm, min_overhang, min_var_reads;
} efilter;

static unsigned char comp_base[256] = {
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32, '!', '"', '#', '$', '%', '&', '\'','(', ')', '*', '+', ',', '-', '.', '/',
  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
  '@', 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
  'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', '[', '\\',']', '^', '_',
  '`', 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
  'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', '{', '|', '}', '~', 127,
  128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
  144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
  176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
  192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
  208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
  224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
  240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

typedef struct {
  int minq;
  double pct;
} read_qual_t;

typedef struct {
  int min_mq, libtype, min_bq, max_depth, all, multi;
  read_qual_t read_qual;
  uint32_t keep_flag[2];
  char *reg, *fai_fname;
  faidx_t *fai;
  regidx_t* idx;
  char **bamfns, **outfns;
  FILE *fps;
  efilter ef;
} mplp_conf_t;

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
  const mplp_conf_t *conf;
  hts_idx_t *idx;
} mplp_aux_t;

typedef struct {
  int n;
  int *n_plp, *m_plp;
  bam_pileup1_t **plp;
} mplp_pileup_t;


SEXP run_sc_pileup(char* bampath,
                   char* fapath,
                   char* region,
                   const char** outfns,
                   int n_outfns,
                   const char** barcodes,
                   int n_bcs,
                   SEXP region_exp,
                   int min_reads,
                   int max_depth,
                   int min_baseQ,
                   int min_mapQ,
                   int libtype,
                   int* b_flags,
                   int* event_filters,
                   double* read_bqual_filter) {

  hts_set_log_level(HTS_LOG_ERROR);

  mplp_aux_t **data;
  int i, tid, *n_plp, ret = 0, n = 1;
  hts_pos_t pos, beg0 = 0, end0 = INT32_MAX, ref_len;

  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  bam_mplp_t iter;
  bam_hdr_t *h = NULL;
  char *ref;

  data = calloc(n, sizeof(mplp_aux_t*));
  plp = calloc(n, sizeof(bam_pileup1_t*));
  n_plp = calloc(n, sizeof(int));

  // read the header of each file in the list and initialize data
  for (i = 0; i < n; ++i) {
    bam_hdr_t *h_tmp;
    data[i] = calloc(1, sizeof(mplp_aux_t));
    data[i]->fp = sam_open(cbampaths[i], "rb");
    if ( !data[i]->fp ) {
      Rf_error("failed to open %s: %s\n", cbampaths[i], strerror(errno));
    }
    if (conf->fai_fname) {
      if (hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
        Rf_error("failed to process %s: %s\n",
                 cfapath, strerror(errno));
      }
    }
    data[i]->conf = conf;
    data[i]->ref = &mp_ref;

    h_tmp = sam_hdr_read(data[i]->fp);
    if ( !h_tmp ) {
      Rf_error("fail to read the header of %s\n", cbampaths[i]);
    }

    if(conf->reg){
      hts_idx_t *idx = NULL;
      idx = sam_index_load(data[i]->fp, cbampaths[i]) ;

      if (idx == NULL) {
        Rf_error("fail to load bamfile index for %s\n", cbampaths[i]);
      }

      if ( (data[i]->iter=sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
        Rf_error("Fail to parse region '%s' with %s\n", conf->reg, cbampaths[i]);
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

  iter = bam_mplp_init(n, readaln, (void**)data);
  // enable overlap detection
  bam_mplp_init_overlaps(iter) ;
  // set max depth
  bam_mplp_set_maxcnt(iter, conf->max_depth);

  if (!iter) {
    REprintf("issue with iterator");
    ret = -1;
    goto fail;
  }

  pcounts *plpc;
  plpc = R_Calloc(n, pcounts);
  for (i = 0; i < n; ++i) {
    // initialize output files
    pd->fps[i] = fopen(R_ExpandFileName(coutfns[i]), "w");
    if (pd->fps[i] == NULL) {
      REprintf("Failed to open file outputfile %s\n", R_ExpandFileName(coutfns[i]));
      ret = -1;
      goto fail;
    }
    plpc[i].pc = R_Calloc(1, counts);
    plpc[i].mc = R_Calloc(1, counts);
  }

  int last_tid = -1;
  int n_iter = 0;
  while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    if (region && (pos < beg0 || pos >= end0)) continue; // not in single region requested

    mplp_get_ref(data[0], tid, &ref, &ref_len);
    if (tid < 0) break;

    // check user interrupt
    if (n_iter % 262144 == 0) {
      if(checkInterrupt()){
        goto fail;
      }
    }
    // ensure position is in requested intervals
    if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, sam_hdr_tid2name(h, tid), pos, pos+1)) continue;

    // get reference base on +/- strand
    int pref_b, mref_b;
    pref_b = (ref && pos < ref_len)? ref[pos] : 'N' ;
    pref_b = toupper(pref_b);
    mref_b = comp_base[(unsigned char) pref_b];

    // reset count structure
    for(i = 0; i < n; ++i){
      clear_pcounts(&plpc[i]);
    }

    // check if site is in a homopolymer
    if(ef->nmer > 0 && check_simple_repeat(&ref, &ref_len, pos, ef->nmer)) continue;

    // check if read count less than min_reads
    int pass_reads = 0;
    for(i = 0; i < n; ++i){
      if (n_plp[i] >= min_reads) {
        pass_reads = 1;
      }
    }
    if(!pass_reads) continue;

    // iterate through bam files
    for (i = 0; i < n; ++i) {
      int j;
      // iterate through reads that overlap position
      for (j = 0; j < n_plp[i]; ++j) {

        const bam_pileup1_t *p = plp[i] + j;

        // get read base
        int c = p->qpos < p->b->core.l_qseq
        ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
        : 'N';

        // store base counts based on library type and r1/r2 status
        int invert = invert_read_orientation(p->b, libtype[i]);
        if(invert < 0){
          REprintf("[internal] invert read orientation failure %i\n", invert);
          ret = -1;
          goto fail;
        }

        // remove bad reads
        int rret = check_read_filters(p, ef, min_baseQ, min_mapQ[i]);
        if(rret > 0) {
          if(rret == 1) {
            if(invert) {
              plpc[i].mc->nx += 1;
            } else {
              plpc[i].pc->nx += 1;
            }
          }
          continue;
        }

        if(invert) c = (char)comp_base[(unsigned char)c];

        // check read for >= mismatch different types and at least n_mm mismatches
        if(ef->n_mm_type > 0 || ef->n_mm > 0) {
          if((invert && mref_b != c) || pref_b != c){
            int m = parse_mismatches(p->b, pos, ef->n_mm_type, ef->n_mm);
            if(m > 0) continue;
          }
        }

        // increment counts
        int cret = count_one_record(p->b, &plpc[i], conf,
                                    sam_hdr_tid2name(h, tid), pos,
                                    pref_b, mref_b,
                                    c, invert, i);
        if(cret < 0) {
          ret = -1;
          goto fail;
        }

      }
    }

    // write or store records if pass depth criteria
    store_counts(pd, plpc, only_keep_variants, n, min_reads, sam_hdr_tid2name(h, tid),
                 pos, pref_b, mref_b,
                 in_mem, ef->min_var_reads);
  }

  fail:
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);

    for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
    free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);

    for (i = 0; i < n; ++i) {
      sam_close(data[i]->fp);
      if(data[i]->iter) hts_itr_destroy(data[i]->iter);
      if(multi_region_itr && data[i]->idx) hts_idx_destroy(data[i]->idx);
      free(data[i]);
      clear_pcounts(&plpc[i]);
      kh_destroy(varhash, plpc[i].mc->var);
      kh_destroy(varhash, plpc[i].pc->var);
      R_Free(plpc[i].pc); R_Free(plpc[i].mc);
      if(!in_mem) fclose(pd->fps[i]);
    }

    free(data); free(plp); free(n_plp);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);

    R_Free(plpc); R_Free(ef);
    if(pd->fps) R_Free(pd->fps);
    R_Free(pd->pdat);
    if (conf->fai) fai_destroy(conf->fai);
    // don't destroy index if passed from R
    if (conf->bed && Rf_isNull(ext)) bed_destroy(conf->bed);

    if (conf->output_reads) {
      clear_rname_set(conf->rnames);
      kh_destroy(rname, conf->rnames);
      fclose(conf->reads_fp);
    }

    if (mismatches) {
      khint_t k;
      for (k = 0; k < kh_end(conf->brhash); ++k){
        if (kh_exist(conf->brhash, k)) free((char*)kh_key(conf->brhash, k));
      }
      kh_destroy(str, conf->brhash);
    }

  if(ret < 0) Rf_error("error detected during pileup");
  return ret;
}

static void check_sc_plp_args(SEXP bampath,
                           SEXP fapath,
                           SEXP region,
                           SEXP region_lst,
                           SEXP event_filters,
                           SEXP min_mapQ,
                           SEXP max_depth,
                           SEXP min_baseQ,
                           SEXP read_bqual_filter,
                           SEXP libtype,
                           SEXP b_flags,
                           SEXP multi_region_itr,
                           SEXP outfns){

  if(!IS_CHARACTER(bampath) || (LENGTH(bampath) != 1)){
    Rf_error("'bampath' must be character(1)");
  }

  if(!IS_CHARACTER(fapath) || (LENGTH(fapath) != 1)){
    Rf_error("'fapath' must be character(1)");
  }

  if(!IS_CHARACTER(region) || (LENGTH(region) > 1)){
    Rf_error("'region' must be character of length 0 or 1");
  }

  if(LENGTH(region_lst) != 3){
    Rf_error("'region_lst' must contain seqnames, pos, and rowidx");
  }

  if(!IS_INTEGER(event_filters) || (LENGTH(event_filters) != 9)){
    Rf_error("'event_filters' must be integer of length 9");
  }

  if(!IS_INTEGER(min_mapQ) || (LENGTH(min_mapQ) != n_files)){
    Rf_error("'min_mapQ' must be integer of same length as bamfiles");
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

  if(!IS_LOGICAL(multi_region_itr) || (LENGTH(multi_region_itr) != 1)){
    Rf_error("'multi_region_itr' must be logical(1)");
  }

  if((!IS_CHARACTER(outfns) || (LENGTH(outfns) != 3))){
    Rf_error("'outfns' must be character(3)");
  }

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

static int set_mplp_conf(mplp_conf_t *conf, int n_bams, char** bamfns,
                         int n_outfns, char** outfns,
                         char* fapath, char* region, int min_mapQ,
                         int min_baseQ, int max_depth, int* b_flags,
                         double* read_bqual_filter, regidx_t* idx,
                         efilter* ef, int* event_filters,
                         int n_bcs, char** bcs){
  int ret = 0;

  int n_bams, char** bamfns
  int n_bams, char** bamfns

  conf->fps = R_Calloc(n_outfns, FILE*);
  if(conf->fps == NULL) Rf_error("[raer internal] unable to alloc");
  conf->outfns = outfns

  if(fapath){
    conf->fai_fname = fapath;
    conf->fai = fai_load(conf->fai_fname);
  }

  if(region) conf->reg = region;
  if(idx) conf->idx = idx;
  conf->min_mq = min_mapQ;

  if(min_baseQ < 0) conf->min_bq = 0;
  if(!max_depth) conf->max_depth = 10000;

  if(b_flags){
    conf->keep_flag[0] = b_flags[0];
    conf->keep_flag[1] = b_flags[1];
  }

  if(read_bqual_filter){
    conf->read_qual.pct = read_bqual_filter[0];
    conf->read_qual.minq = (int)read_bqual_filter[1];
  }

  if(nbcs > 0){
    void *cbhash = khash_str2int_init();
    for(int i = 0; i < nbcs; ++i){
      char *cb = strdup(bcs[i]);
      if(khash_str2int_inc(cbhash, cb) == -1){
        Rf_error("[raer internal] unable to populate barcode hashmap");
      }
    }
    conf->cbhash = cbhash;
  } else {
    conf->cbhash = NULL;
  }

  conf->ef = set_event_filters(ef, event_filters);
  return(ret);
}
SEXP do_run_pileup(SEXP bampath,
                   SEXP fapath,
                   SEXP region,
                   SEXP region_lst,
                   SEXP region_exp,
                   SEXP event_filters,
                   SEXP min_mapQ,
                   SEXP max_depth,
                   SEXP min_baseQ,
                   SEXP read_bqual_filter,
                   SEXP libtype,
                   SEXP b_flags,
                   SEXP outfns) {

  check_sc_plp_args(bampath, fapath, region,
                    region_lst,
                    event_filters,
                    min_mapQ,
                    max_depth,
                    min_baseQ,
                    read_bqual_filter,
                    libtype,
                    b_flags,
                    multi_region_itr,
                    outfns);

  int i;
  char * cbampath = (char *) translateChar(STRING_ELT(bampath, 0));
  char * cfapath = (char *) translateChar(STRING_ELT(fapath, 0));
  char * cregion = LENGTH(region) == 0 ?
                   NULL : (char *) translateChar(STRING_ELT(region, 0));

  const char** coutfns;
  coutfns = (const char**) R_alloc(sizeof(const char *), Rf_length(outfns));
  for (i = 0; i < LENGTH(outfns); ++i){
    coutfns[i] = (char *) translateChar(STRING_ELT(outfns, i));
  }

  mplp_conf_t mplp, *conf;
  memset(&mplp, 0, sizeof(mplp_conf_t));
  conf = &mplp;

  if(!Rf_isNull(region_exp)){
    _REG_IDX *ridx = REGIDX(region_exp) ;
    if (ridx->index == NULL){
      Rf_error("Failed to load region index");
    }
    conf->idx =ridx->index;
  } else if(!Rf_isNull(region_lst)) {
    SEXP region_exp = regidx_build(region_lst);
    REG_IDX *ridx = REGIDX(region_exp) ;
    if (ridx->index == NULL){
      Rf_error("Failed to load region index");
    }
    conf->idx = ridx->index;
  } else {
    conf->idx = NULL;
  }

  int ret, res;
  efilter* ef = R_Calloc(1, efilter);
  ret = set_mplp_conf(conf, ridx->index, ef, INTEGER(event_filters));

  res = run_cpileup(cbampaths,);                                                                                                            ext);

  return res ;
}
