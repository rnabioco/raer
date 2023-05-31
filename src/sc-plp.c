#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/hts_log.h>
#include <htslib/regidx.h>
#include <htslib/khash_str2int.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "regfile.h"
#include "plp_utils.h"

typedef struct umi_bases_t {
    int ref;
    int alt;
    int oth;
} umi_bases_t ;

// hashmap with UMI seq as key pointing to a struct storing base counts
KHASH_MAP_INIT_STR(umimap, umi_bases_t*)
typedef khash_t(umimap)* umimap_t;

typedef struct  {
  umimap_t umi; // hashmap storing umi sequences, and base count array
  int ref; // # of UMIs supporting ref
  int alt; // # of UMIs supporting alt
} cb_t;

KHASH_MAP_INIT_STR(cbumimap, cb_t*)
typedef khash_t(cbumimap)* cbumi_map_t;


static void clear_umi(umimap_t uhash) {
  khint_t k;
  umi_bases_t* ub;

  if (!uhash) return;
  for (k = kh_begin(uhash); k < kh_end(uhash); ++k) {
    if (kh_exist(uhash, k)) {
        free((char*)kh_key(uhash, k));
        ub = (umi_bases_t*) kh_val(uhash, k);
        if(ub) free(ub);
    }
  }
  kh_clear(umimap, uhash);
}

static void clear_cb_umiset(cbumi_map_t cbhash) {
  khint_t k;
  cb_t* cdat;
  if (!cbhash) return;
  for (k = kh_begin(cbhash); k < kh_end(cbhash); ++k) {
    if (!kh_exist(cbhash, k)) continue;
    cdat = kh_value(cbhash, k);
    clear_umi(cdat->umi);
    cdat->ref = cdat->alt = 0;
  }
}

static cb_t* init_umihash() {
  cb_t* cb = R_Calloc(1, cb_t);
  cb->umi = kh_init(umimap);
  return cb ;
}

static void free_hashmaps(cbumi_map_t cbhash, str2intmap_t cbidx) {
  khint_t k;
  cb_t* cdat;
  if (cbhash) {
    for (k = kh_begin(cbhash); k < kh_end(cbhash); ++k) {
      if (!kh_exist(cbhash, k)) continue;
      cdat = kh_value(cbhash, k);
      clear_umi(cdat->umi);
      kh_destroy(umimap, cdat->umi);
      free(cdat);
    }
    kh_destroy(cbumimap, cbhash);
  }

  if (cbidx) {
    for (k = kh_begin(cbidx); k < kh_end(cbidx); ++k) {
      if (kh_exist(cbidx, k)) free((char*)kh_key(cbidx, k));
    }
    kh_destroy(str2intmap, cbidx);
  }

}

// single cell pileup parameters
typedef struct {
  int min_mq, libtype, min_bq, max_depth;
  read_qual_t read_qual;
  uint32_t keep_flag[2];
  char* qregion;
  regidx_t* reg_idx;
  regitr_t* reg_itr;
  char* mtxfn;
  char* bcfn;
  char* sitesfn;
  FILE** fps; // [0] = mtxfn, [1] = sitefn, [2] = bcfn;
  efilter ef;
  str2intmap_t cbidx; // hashmap cellbarcode key -> column index value
  cbumi_map_t cbmap; // hashmap cellbarcode key -> cbumi_map_t as value
  int has_umi;
  char* umi_tag;
  char* cb_tag;
  int idx_skip;
  int pe; // 1 if paired end, 0 if not
  int min_counts; // if 0 report all sites in sparseMatrix
  int site_idx; // counter of sites written, used for row index if min_counts > 0
  int is_ss2; // 1 if smart-seq2 mode, 0 otherwise
} sc_mplp_conf_t;

typedef struct {
  samFile* fp;
  hts_itr_t* iter;
  bam_hdr_t* h;
  const sc_mplp_conf_t* conf;
  hts_idx_t* idx;
} mplp_sc_aux_t;


/* return codes,
 * 0 = read passes,
 * 1 = read fails
 * 2 = read fails due to refskip, deletion, overlapping mate
 */
static int check_read_filters(const bam_pileup1_t* p, sc_mplp_conf_t* conf) {

  // skip indel and ref skip ;
  if (p->is_del || p->is_refskip) return (2) ;

  // filter based on mapq
  if (p->b->core.qual < conf->min_mq) return (1) ;

  // check base quality
  int bq = p->qpos < p->b->core.l_qseq
           ? bam_get_qual(p->b)[p->qpos]
           : 0;
  if (bq < conf->min_bq) {
    if (bq == 0) { // overlapping read pairs get base qualities set to 0 by mplp
      return (2);
    } else {
      return (1);
    }
  }

  // check if pos is within x dist from 5' end of read, qpos is 0-based
  if (check_variant_pos(p->b, p->qpos, conf->ef.trim_5p_dist, conf->ef.trim_3p_dist)) return (1);

  // check for splice in alignment nearby
  if (conf->ef.splice_dist && dist_to_splice(p->b, p->qpos, conf->ef.splice_dist) >= 0) return (1);

  // check if site in splice overhang and > min_overhang
  if (conf->ef.min_overhang && check_splice_overhang(p->b, p->qpos, conf->ef.min_overhang) > 0) return (1);

  // check if indel event nearby
  if (conf->ef.indel_dist && dist_to_indel(p->b, p->qpos, conf->ef.indel_dist) >= 0) return (1);

  return 0;
}

/* return codes,
 * -1 error in hashmap
 * 0 missing cb or umi
 * 1 umi is duplicate, or base/strand not in payload
 * 2 umi is reference
 * 3 umi is alternate
 */
static int count_record(bam1_t* b, sc_mplp_conf_t* conf, payload_t* pld,
                        unsigned char base, int strand, char* bamid) {
  khiter_t k;
  char* cb, *cb_cpy, *umi, *umi_val;
  cb_t* cbdat;
  umi_bases_t* ucounts;

  int cret, uret;
  if (conf->is_ss2) {
    cb = bamid;
  } else {
    cb = get_aux_ztag(b, conf->cb_tag);
    if (cb == NULL) return (0);
  }

  if (pld->strand != strand) return (1);

  cb_cpy = strdup(cb);

  k = kh_get(str2intmap, conf->cbidx, cb_cpy);
  if (k == kh_end(conf->cbidx)) {
    free(cb_cpy);
    return (0);
  }

  k = kh_put(cbumimap, conf->cbmap, cb_cpy, &cret);
  if (cret < 0) {
    free(cb_cpy);
    return (-1);
  }
  if (cret == 0) free(cb_cpy); // CB in hash already
  if (cret == 1) kh_value(conf->cbmap, k) = init_umihash();
  cbdat = kh_value(conf->cbmap, k);

  if (conf->has_umi) {
    umi = get_aux_ztag(b, conf->umi_tag);
    if (umi == NULL) return (0);
    umi_val = strdup(umi);
    k = kh_put(umimap, cbdat->umi, umi_val, &uret);

    if (uret < 0) {
      free(umi_val);
      return (-1);
    }
    if (uret == 0) free(umi_val); // umi in hash already
    if (uret == 1) kh_value(cbdat->umi, k) = R_Calloc(1, umi_bases_t);

    ucounts = kh_value(cbdat->umi, k);

    if ((unsigned char)*pld->ref == base) {
        ucounts->ref++;
        return (2);
    } else if ((unsigned char)*pld->alt == base) {
        ucounts->alt++;
        return (3);
    } else {
        ucounts->oth++;
        return(1);
    }
  } else {
      if ((unsigned char)*pld->ref == base) {
          cbdat->ref++;
          return (2);
      } else if ((unsigned char)*pld->alt == base) {
          cbdat->alt++;
          return (3);
      }
  }
  return (1);
}

static int count_consensus_bases(cb_t* cbmap) {
  khint_t k;
  umi_bases_t* uc;

  for (k = kh_begin(cbmap->umi); k != kh_end(cbmap->umi); ++k) {
    if (!kh_exist(cbmap->umi, k)) continue;
    uc = kh_val(cbmap->umi, k);
    if(!uc) continue;

    // umi supports other bases
    if((uc->oth > uc->alt) && (uc->oth > uc->ref)) continue;

    // probably shouldn't happen
    if((uc->alt == 0) && (uc->ref == 0)) continue;

    // in case of ties, report as alt
    if(uc->alt >= uc->ref) {
        cbmap->alt++;
    } else {
        cbmap->ref++;
    }
  }

  return(0);
}

/* write entry to sparse array-like format
 * col 1 = row index (site index)
 * col 2 = column index (barcode index)
 * col 3 = values for ref
 * col 4 = values for alt
 *
 * read into R as a 4 column matrix, then coerce to list of 2 sparseMatrices.
 */
static int write_counts(sc_mplp_conf_t* conf, payload_t* pld, const char* seqname, int pos) {
  const char* cb;
  cb_t* cbdat;
  int c_idx, r_idx, ret, n_rec = 0;
  khint_t k, j;

  for (k = kh_begin(conf->cbmap); k != kh_end(conf->cbmap); ++k) {
    if (!kh_exist(conf->cbmap, k)) continue;
    cb = kh_key(conf->cbmap, k);
    cbdat = kh_val(conf->cbmap, k);

    j = kh_get(str2intmap, conf->cbidx, cb);
    if (j != kh_end(conf->cbidx)) {
      c_idx = kh_value(conf->cbidx, j);
    } else {
      REprintf("[raer internal] error retrieving CB %s %d\n", cb, j);
      return (-1);
    }

    if (conf->has_umi) {
      if(!cbdat->umi) {
          REprintf("[raer internal] no bases found for CB %s %d\n", cb, j);
          return (-1);
      }
      if(count_consensus_bases(cbdat) < 0) {
          REprintf("[raer internal] error calculating consensus base for umi\n");
          return (-1);
      }
    }
    if (cbdat->ref == 0 && cbdat->alt == 0) continue;
    r_idx = conf->min_counts == 0 ? pld->idx : conf->site_idx;
    ret = fprintf(conf->fps[0], "%d %d %d %d\n", r_idx, c_idx, cbdat->ref, cbdat->alt);
    if (ret < 0) return (-1);
    n_rec += 1;
  }

  // write site
  if (conf->min_counts > 0) {
    ret = fprintf(conf->fps[1], "%d\t%s\t%d\t%d\t%s\t%s\n", pld->idx, seqname, pos + 1, pld->strand, pld->ref, pld->alt);
    if (ret < 0) return (-1);
    conf->site_idx += 1;
  }
  return (n_rec);
}

// read processing function for pileup
static int screadaln(void* data, bam1_t* b) {
  mplp_sc_aux_t* g = (mplp_sc_aux_t*)data;
  int ret, skip = 0;
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
    if (g->conf->cbidx && g->conf->cb_tag) {
      cb = get_aux_ztag(b, g->conf->cb_tag);
      if (!cb) {
        skip = 1;
        continue;
      }
      k = kh_get(str2intmap, g->conf->cbidx, cb);
      if (k == kh_end(g->conf->cbidx)) {
        skip = 1;
        continue;
      }
    }

    // test required and filter flags
    test_flag = (g->conf->keep_flag[0] & ~b->core.flag) |
                (g->conf->keep_flag[1] & b->core.flag);
    if (~test_flag & 2047u) {
      skip = 1;
      continue;
    }

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
    } else if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) {
      skip = 1;
      // check if read quality is > x in at least y% of read
    } else if (g->conf->read_qual.pct &&
               g->conf->read_qual.minq &&
               !read_base_quality(b, g->conf->read_qual.pct, g->conf->read_qual.minq)) {
      skip = 1;
    }
  } while (skip);
  return ret;
}

static void set_event_filters(efilter* ef, int* event_filters) {
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

// process one site
static int process_one_site(const bam_pileup1_t** plp, sc_mplp_conf_t* conf,
                            payload_t* pld, bam_hdr_t* h, int tid,
                            int pos, int nbam, int* n_plp, char* bamid) {
    int i, j,n_ref,n_alt;
    i = j = n_ref = n_alt = 0;

    clear_cb_umiset(conf->cbmap);

    for (i = 0; i < nbam; ++i) {
        // iterate through reads that overlap position
        for (j = 0; j < n_plp[i]; ++j) {

            const bam_pileup1_t* p = plp[i] + j;

            // get read base
            int c = p->qpos < p->b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
            : 'N';

            // determine if based should be complemented based on library type and r1/r2 status
            int invert = invert_read_orientation(p->b, conf->libtype);
            if (invert < 0) {
                REprintf("[raer internal] invert read orientation failure %i\n", invert);
                return(-1);
            }

            // remove bad reads
            int rret = check_read_filters(p, conf);
            if (rret > 0) continue;

            if (invert) c = (char)comp_base[(unsigned char)c];

            // invert is 0 if pos, 1 if neg
            // set strand to match payload with 1 is pos 2 is neg
            int strand = invert + 1; //

            int cret = count_record(p->b, conf, pld, c, strand, bamid);

            if(cret < 0) {
                REprintf("[raer internal] issue with counting records\n");
                return(-1);
            } else if (cret == 2) {
                n_ref += 1;
            } else if (cret == 3) {
                n_alt += 1;
            }
        }
    }

    // write records
    int n_counted = n_ref + n_alt;
    if (conf->min_counts == 0 ||
        (n_counted >= conf->min_counts && n_alt >= conf->ef.min_var_reads) ) {

        int wr = write_counts(conf, pld, sam_hdr_tid2name(h, tid), pos);
        if (wr == -1) {
            REprintf("[raer internal] error in writing output.");
            return(-1);
        }
    }
    return(1);
}

static int run_scpileup(sc_mplp_conf_t* conf, char* bamfn, char* index, char* bamid) {

  hts_set_log_level(HTS_LOG_ERROR);

  mplp_sc_aux_t** data;
  int i, tid, *n_plp, ret = 0;
  hts_pos_t pos, beg0 = 0, end0 = INT32_MAX;

  const bam_pileup1_t** plp;
  bam_mplp_t iter = NULL;
  bam_hdr_t* h = NULL;

  int nbam  = 1;
  data = calloc(1, sizeof(mplp_sc_aux_t*));
  plp = calloc(1, sizeof(bam_pileup1_t*));
  n_plp = calloc(1, sizeof(int));

  // read the header of bam file and initialize data

  data[0] = calloc(1, sizeof(mplp_sc_aux_t));
  data[0]->fp = sam_open(bamfn, "rb");
  if (!data[0]->fp) {
    REprintf("[raer interal] failed to open %s: %s\n", bamfn, strerror(errno));
    ret = -1;
    goto fail;
  }

  data[0]->conf = conf;

  h = sam_hdr_read(data[0]->fp);
  if (!h) {
    REprintf("[raer interal] fail to read the header of %s\n", bamfn);
    ret = -1;
    goto fail;
  }

  if (conf->qregion) {
    hts_idx_t* idx = NULL;
    idx = sam_index_load2(data[0]->fp, bamfn, index) ;
    if (idx == NULL) {
      REprintf("[raer interal] fail to load bamfile index for %s\n", bamfn);
      ret = -1;
      goto fail;
    }

    if ((data[0]->iter=sam_itr_querys(idx, h, conf->qregion)) == 0) {
      REprintf("[raer interal] Fail to parse region '%s' with %s\n", conf->qregion, bamfn);
      ret = -1;
      goto fail;
    }

    beg0 = data[0]->iter->beg;
    end0 = data[0]->iter->end;
    hts_idx_destroy(idx);

  } else {
    data[0]->iter = NULL;
  }

  data[0]->h = h;

  iter = bam_mplp_init(nbam, screadaln, (void**)data);
  // enable overlap detection
  if (conf->pe) bam_mplp_init_overlaps(iter) ;
  // set max depth
  bam_mplp_set_maxcnt(iter, conf->max_depth);

  if (!iter) {
    REprintf("[raer interal] issue with iterator");
    ret = -1;
    goto fail;
  }

  for (i = 0; i < 3; ++i) {
    if (conf->fps[i] == NULL) {
      REprintf("[raer interal] Failed to open file outputfiles\n");
      ret = -1;
      goto fail;
    }
  }

  int n_iter = 0;
  while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    if (conf->qregion && (pos < beg0 || pos >= end0)) continue; // not in single region requested
    if (tid < 0) break;

    // check user interrupt
    if (n_iter % 65536 == 0) {
      if (checkInterrupt()) goto fail;
    }

    // ensure position is in requested intervals
    // if so recover payload with information about strand and ref/alt bases
    int ol;
    payload_t* pld;
    if (!(conf->reg_idx && tid >= 0)) continue;

    ol = regidx_overlap(conf->reg_idx,
                          sam_hdr_tid2name(h, tid),
                          pos,
                          pos,
                          conf->reg_itr);
    if (!ol) continue;

    // iterate through multiple requested alleles per position/region
    while ( regitr_overlap(conf->reg_itr) ) {
      pld = regitr_payload(conf->reg_itr, payload_t*);
      ret = process_one_site(plp, conf, pld, h, tid, pos, nbam, n_plp, bamid);
      if(ret < 0) goto fail;
      ++n_iter;
    }
  }

fail:
  if (iter) bam_mplp_destroy(iter);
  bam_hdr_destroy(h);

  for (i = 0; i < nbam; ++i) {
    sam_close(data[0]->fp);
    if (data[0]->iter) hts_itr_destroy(data[0]->iter);
    if (data[0]->idx) hts_idx_destroy(data[0]->idx);
    free(data[0]);
  }
  free(data);
  free(plp);
  free(n_plp);

  return ret;
}


static int set_sc_mplp_conf(sc_mplp_conf_t* conf, int nbams,
                            int n_outfns, char** outfns, char* qregion, regidx_t* idx,
                            int min_mapQ, int min_baseQ, double* read_bqual_filter,
                            int max_depth, int* b_flags, int* event_filters, int libtype,
                            int n_bcs, char** bcs, char* cbtag, char* umi,
                            int idx_skip, int pe, int min_counts) {
  conf->is_ss2 = nbams > 1 ? 1 : 0;

  conf->fps = R_Calloc(n_outfns, FILE*);
  if (conf->fps == NULL) {
    REprintf("[raer internal] unable to alloc");
    return (-1);
  }
  conf->mtxfn = outfns[0];
  conf->sitesfn = outfns[1];
  conf->bcfn = outfns[2];

  // initialize output files
  conf->fps[0] = fopen(R_ExpandFileName(conf->mtxfn), "w");
  conf->fps[1] = fopen(R_ExpandFileName(conf->sitesfn), "w");
  conf->fps[2] = fopen(R_ExpandFileName(conf->bcfn), "w");

  if (qregion) conf->qregion = qregion;
  if (idx) {
    conf->reg_idx = (regidx_t*)idx;
    conf->reg_itr = regitr_init(conf->reg_idx);
  }

  conf->min_mq = (min_mapQ < 0) ? 0 : min_mapQ;
  conf->min_bq = (min_baseQ < 0) ? 0 : min_baseQ;
  conf->max_depth = (!max_depth) ? 10000: max_depth;

  if (b_flags) {
    conf->keep_flag[0] = b_flags[0];
    conf->keep_flag[1] = b_flags[1];
  }

  if (read_bqual_filter) {
    conf->read_qual.pct = read_bqual_filter[0];
    conf->read_qual.minq = (int)read_bqual_filter[1];
  }

  memset(&conf->ef, 0, sizeof(efilter));
  set_event_filters(&(conf->ef), event_filters);

  conf->libtype = libtype;
  int hret = 0;
  khint_t k;
  if (n_bcs > 0) {
    conf->cbmap = kh_init(cbumimap);
    conf->cbidx = kh_init(str2intmap);

    for (int i = 0; i < n_bcs; ++i) {
      char* cb = strdup(bcs[i]);
      k = kh_put(str2intmap, conf->cbidx, cb, &hret);
      if (hret == -1 || hret == 2) {
        REprintf("[raer internal] unable to populate barcode hashmap");
        return (-1);
      } else if (hret == 0) {
        free(cb); // was a duplicate cellbarcode
      } else {
        // storing cb index as 1-based (for compatibility with matrix-market format)
        kh_value(conf->cbidx, k) = i + 1;
      }
    }
  } else {
    conf->cbidx = NULL;
  }

  if (umi) {
    conf->has_umi = 1;
    conf->umi_tag = umi;
  } else {
    conf->has_umi = 0;
    conf->umi_tag = NULL;
  }
  conf->cb_tag = cbtag;
  conf->idx_skip = idx_skip;
  conf->pe = pe;
  conf->min_counts = min_counts < 0 ? 0 : min_counts;
  conf->site_idx = 1;

  return (0);
}

static void write_barcodes(FILE* fp, char** bcs, int n) {
  for (int i = 0; i < n; ++i) {
    fprintf(fp, "%s\n", bcs[i]);
  }
}

static int write_all_sites(sc_mplp_conf_t* conf) {
  if (!conf->reg_idx) return (-1);
  regitr_t* itr = regitr_init(conf->reg_idx);
  payload_t* pld;
  while (regitr_loop(itr)) {
    pld = regitr_payload(itr, payload_t*);
    fprintf(conf->fps[1], "%d\t%s\t%d\t%d\t%s\t%s\n",
            pld->idx,
            itr->seq,
            (int)itr->beg+1,
            pld->strand,
            pld->ref,
            pld->alt);

  }
  return (1);
}


static void check_sc_plp_args(SEXP bampaths, SEXP indexes, SEXP qregion, SEXP lst,
                              SEXP barcodes, SEXP cbtag, SEXP event_filters,
                              SEXP min_mapQ, SEXP max_depth, SEXP min_baseQ,
                              SEXP read_bqual_filter, SEXP libtype, SEXP b_flags,
                              SEXP outfns, SEXP umi, SEXP index_skip, SEXP pe,
                              SEXP min_counts) {

  if (!IS_CHARACTER(bampaths) || (LENGTH(bampaths) < 1)) {
    Rf_error("'bampaths' must be character");
  }

  if (!IS_CHARACTER(indexes) || (LENGTH(indexes) < 1)) {
      Rf_error("'indexes' must be character");
  }

  if (LENGTH(bampaths) != LENGTH(indexes)) {
      Rf_error("'bampaths' and 'indexes' must be same length");
  }

  if (!IS_CHARACTER(qregion) || (LENGTH(qregion) > 1)) {
    Rf_error("'qregion' must be character of length 0 or 1");
  }

  if (Rf_isNull(lst)) {
    Rf_error("'lst' must be non-null");
  }

  if (!IS_CHARACTER(barcodes) || (LENGTH(barcodes) < 1)) {
    Rf_error("'barcodes' must be character");
  }

  if ((!IS_CHARACTER(cbtag) || (LENGTH(cbtag) > 1))) {
    Rf_error("'cbtag' must be character of length 0 or 1");
  }

  if (!IS_INTEGER(event_filters) || (LENGTH(event_filters) != 9)) {
    Rf_error("'event_filters' must be integer of length 9");
  }

  if (!IS_INTEGER(min_mapQ) || (LENGTH(min_mapQ) != 1)) {
    Rf_error("'min_mapQ' must be integer(1)");
  }

  if (!IS_INTEGER(max_depth) || (LENGTH(max_depth) != 1)) {
    Rf_error("'max_depth' must be integer(1)");
  }

  if (!IS_INTEGER(min_baseQ) || (LENGTH(min_baseQ) != 1)) {
    Rf_error("'min_baseQ' must be integer(1)");
  }

  if (!IS_NUMERIC(read_bqual_filter) || (LENGTH(read_bqual_filter) != 2)) {
    Rf_error("'read_bqual_filter' must be numeric of length 2");
  }

  if (!IS_INTEGER(libtype) || (LENGTH(libtype) != 1)) {
    Rf_error("'lib_type' must be integer(1)");
  }

  if (!IS_INTEGER(b_flags) || (LENGTH(b_flags) != 2)) {
    Rf_error("'b_flags' must be integer of length 2");
  }

  if ((!IS_CHARACTER(outfns) || (LENGTH(outfns) != 3))) {
    Rf_error("'outfns' must be character(3)");
  }

  if ((!IS_CHARACTER(umi) || (LENGTH(umi) > 1))) {
    Rf_error("'umi' must be character of length 0 or 1");
  }

  if (!IS_LOGICAL(index_skip) || (LENGTH(index_skip) != 1)) {
    Rf_error("'index_skip' must be logical(1)");
  }

  if (!IS_LOGICAL(pe) || (LENGTH(pe) != 1)) {
    Rf_error("'pe' must be logical(1)");
  }

  if (!IS_INTEGER(min_counts) || (LENGTH(min_counts) != 1)) {
    Rf_error("'min_counts' must be integer(1)");
  }

}

SEXP scpileup(SEXP bampaths, SEXP indexes, SEXP query_region, SEXP lst,
              SEXP barcodes, SEXP cbtag, SEXP event_filters, SEXP min_mapQ,
              SEXP max_depth, SEXP min_baseQ, SEXP read_bqual_filter,
              SEXP libtype, SEXP b_flags, SEXP outfns, SEXP umi,
              SEXP index_skip, SEXP pe, SEXP min_counts) {

  check_sc_plp_args(bampaths, indexes, query_region, lst,
                    barcodes, cbtag, event_filters, min_mapQ, max_depth,
                    min_baseQ, read_bqual_filter, libtype,
                    b_flags, outfns, umi, index_skip, pe, min_counts);

  regidx_t* idx = regidx_build(lst, 1);
  if (!idx) Rf_error("Failed to build region index");

  int i, ret = 0, res = 0;
  char** cbampaths;
  char** cindexes;
  int nbams = Rf_length(bampaths);
  cbampaths = (char**) R_alloc(sizeof(char*), nbams);
  cindexes = (char**) R_alloc(sizeof(char*), nbams);
  for (i = 0; i < nbams; ++i) {
    cbampaths[i] = (char*) translateChar(STRING_ELT(bampaths, i));
    cindexes[i]  = (char*) translateChar(STRING_ELT(indexes, i));
  }

  char* cq_region = LENGTH(query_region) == 0 ?
                    NULL : (char*) translateChar(STRING_ELT(query_region, 0));

  char** coutfns;
  int nout = Rf_length(outfns);
  coutfns = (char**) R_alloc(sizeof(char*), nout);
  for (i = 0; i < nout; ++i) {
    coutfns[i] = (char*) translateChar(STRING_ELT(outfns, i));
  }

  char** bcs;
  int nbcs = Rf_length(barcodes);
  bcs = (char**) R_alloc(sizeof(char*), nbcs);
  for (i = 0; i < nbcs; ++i) {
    bcs[i] = (char*) translateChar(STRING_ELT(barcodes, i));
  }

  char* c_cbtag = LENGTH(cbtag) == 0 ?
                  NULL : (char*) translateChar(STRING_ELT(cbtag, 0));

  char* c_umi = LENGTH(umi) == 0 ?
                NULL : (char*) translateChar(STRING_ELT(umi, 0));

  sc_mplp_conf_t ga;
  memset(&ga, 0, sizeof(sc_mplp_conf_t));

  ret = set_sc_mplp_conf(&ga, nbams, nout, coutfns,
                         cq_region, idx, INTEGER(min_mapQ)[0],
                         INTEGER(min_baseQ)[0], REAL(read_bqual_filter),
                         INTEGER(max_depth)[0], INTEGER(b_flags),
                         INTEGER(event_filters), INTEGER(libtype)[0],
                         nbcs, bcs, c_cbtag, c_umi, LOGICAL(index_skip)[0],
                         LOGICAL(pe)[0], INTEGER(min_counts)[0]);

  if (ret >= 0) {
    // write barcodes file, all barcodes will be reported in matrix
    write_barcodes(ga.fps[2], bcs, nbcs);
    // write sites, if all sites requested, otherwise write during pileup
    if (ga.min_counts == 0) write_all_sites(&ga);
    if (ga.is_ss2) {
      for (i = 0; i < nbams; ++i) {
        res = run_scpileup(&ga, cbampaths[i], cindexes[i], bcs[i]);
        if (res < 0) {
          REprintf("[raer internal] error processing bamfile %s:", cbampaths[i]);
          break;
        }
      }
    } else {
      res = run_scpileup(&ga, cbampaths[0], cindexes[0], NULL);
    }
  }

  for (i = 0; i < nout; ++i) {
    if (ga.fps[i]) fclose(ga.fps[i]);
  }

  if (ga.fps) R_Free(ga.fps);
  if (ga.reg_itr) regitr_destroy(ga.reg_itr);
  if (ga.reg_idx) regidx_destroy(ga.reg_idx);
  free_hashmaps(ga.cbmap, ga.cbidx);

  if (res < 0) REprintf("[raer internal] error detected during pileup, %d\n", res);

  return ScalarInteger(res) ;
}


