#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "bedidx.h"
#include "bedfile.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <inttypes.h>
#include <Rinternals.h>

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
  int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, all, rev_del;
  int rflag_require, rflag_filter;
  char *reg, *pl_list, *fai_fname, *output_fname;
  faidx_t *fai;
  void *bed, *rghash, *auxlist;
  int argc;
  char **argv;
  char sep, empty, no_ins, no_ins_mods, no_del, no_ends;
} mplp_conf_t;

typedef struct {
  char *ref[2];
  int ref_id[2];
  int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
  samFile *fp;
  hts_itr_t *iter;
  bam_hdr_t *h;
  mplp_ref_t *ref;
  const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
  int n;
  int *n_plp, *m_plp;
  bam_pileup1_t **plp;
} mplp_pileup_t;

static int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, int *ref_len) {
  mplp_ref_t *r = ma->ref;

  //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

  if (!r || !ma->conf->fai) {
    *ref = NULL;
    return 0;
  }

  // Do we need to reference count this so multiple mplp_aux_t can
  // track which references are in use?
  // For now we just cache the last two. Sufficient?
  if (tid == r->ref_id[0]) {
    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
  }
  if (tid == r->ref_id[1]) {
    // Last, swap over
    int tmp_id;
    int tmp_len;
    tmp_id  = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp_id;
    tmp_len = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp_len;

    char *tc;
    tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
  }

  // New, so migrate to old and load new
  free(r->ref[1]);
  r->ref[1]     = r->ref[0];
  r->ref_id[1]  = r->ref_id[0];
  r->ref_len[1] = r->ref_len[0];

  r->ref_id[0] = tid;
  r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                               ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

  if (!r->ref[0]) {
    r->ref[0] = NULL;
    r->ref_id[0] = -1;
    r->ref_len[0] = 0;
    *ref = NULL;
    return 0;
  }

  *ref = r->ref[0];
  *ref_len = r->ref_len[0];
  return 1;
}

// read level processing for pileup
static int readaln(void *data, bam1_t *b) {
  mplp_aux_t *g = (mplp_aux_t *)data;
  int ret;
  int skip = 0;
  while (1) {
    ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->h, b);
    if (ret < 0) break;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

    if (g->conf->bed) { // test overlap
      skip = !bed_overlap(g->conf->bed, g->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
      if (skip) continue;
    }
    // check for mapping status
    uint8_t* aux_info = bam_aux_get(b, "NH") ;
    int n_errors = 0 ;
    if(!aux_info) {
      if(n_errors < 5){
        REprintf("Warning: missing tag info in reads, e.g: %s", bam_get_qname(b)) ;
        n_errors += 1 ;
      }
      continue ;
    }

    int64_t mmap_tag = bam_aux2i(aux_info);
    if(mmap_tag > 1) continue ;

    // check for proper pair, if paired end
    if((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) continue ;
    break;
  }

  return ret;
}

// struct for holding counts
typedef struct  {
  int pos;
  int ptotal, pnr, pnv, pna, pnt, png, pnc, pnn;
  int mtotal, mnr, mnv, mna, mnt, mng, mnc, mnn;
  int pref_b;
  int mref_b;
} pcounts;


static void print_plus_counts(FILE *fp, pcounts *pc, int n, char* ctig, int pos, int pref){
  fprintf(fp,
          "%s\t%i\t%c\t%c",
          ctig,
          pos + 1,
          '+',
          pref) ;
  int i;
  for(i = 0; i < n; ++i){
    fprintf(fp,
            "\t%i\t%i\t%i\t%i\t%i\t%i\t%i",
            pc[i].pnr,
            pc[i].pnv,
            pc[i].pna,
            pc[i].pnt,
            pc[i].pnc,
            pc[i].png,
            pc[i].pnn);
  }
  fprintf(fp,"\n");
}

static void print_minus_counts(FILE *fp, pcounts *pc, int n, char* ctig, int pos, int mref){
  fprintf(fp,
          "%s\t%i\t%c\t%c",
          ctig,
          pos + 1,
          '-',
          mref) ;
  int i;
  for(i = 0; i < n; ++i){
    fprintf(fp,
            "\t%i\t%i\t%i\t%i\t%i\t%i\t%i",
            pc[i].mnr,
            pc[i].mnv,
            pc[i].mna,
            pc[i].mnt,
            pc[i].mnc,
            pc[i].mng,
            pc[i].mnn);
  }
  fprintf(fp,"\n");
}

static void print_counts(FILE *fp, pcounts *pc, int n, int min_1, int min_2,
                         char* ctig, int pos, int pref, int mref){

  // determine if counts pass thresholds
  if(n == 1){
    if(pc[0].ptotal >= min_1){
      print_plus_counts(fp, pc, n, ctig, pos, pref) ;
    }

    if(pc[0].mtotal >= min_1){
      print_minus_counts(fp, pc, n, ctig, pos, mref) ;
    }
  } else if (n == 2) {

    if(pc[0].ptotal >= min_1 && pc[1].ptotal >= min_2){
      print_plus_counts(fp, pc, n, ctig, pos, pref) ;
    }

    if(pc[0].mtotal >= min_1 && pc[1].mtotal >= min_2){
      print_minus_counts(fp, pc, n, ctig, pos, mref) ;
    }
  }
}


int run_cpileup(const char** cbampaths,
                int n,
                char* cfapath,
                char* cregion,
                char* coutfn,
                char* cbedfn,
                int* min_reads,
                int max_depth,
                int min_baseQ,
                int libtype,
                SEXP ext) {

  if (n > 2 || n < 1) {
    REprintf("pileup requires 1 or 2 bam files");
    return 1;
  }

  mplp_aux_t **data;
  int i, tid, *n_plp, tid0 = 0;
  int pos, beg0 = 0, end0 = INT32_MAX, ref_len;

  if(!min_baseQ) min_baseQ = 20;
  if(!max_depth) max_depth = 10000;

  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  bam_mplp_t iter;
  bam_hdr_t *h = NULL; /* header of first file in input list */
  char *ref;
  FILE *pileup_fp = NULL;

  mplp_pileup_t gplp;

  mplp_conf_t mplp;
  mplp_conf_t *conf;
  memset(&mplp, 0, sizeof(mplp_conf_t));
  conf = &mplp;

  memset(&gplp, 0, sizeof(mplp_pileup_t));
  data = calloc(n, sizeof(mplp_aux_t*));
  plp = calloc(n, sizeof(bam_pileup1_t*));
  n_plp = calloc(n, sizeof(int));

  if(cfapath){
    conf->fai_fname = cfapath;
    conf->fai = fai_load(conf->fai_fname);
  }

  if(cbedfn){
    conf->bed = bed_read(cbedfn);
  } else if (!Rf_isNull(ext)){
    _BED_FILE *ffile = BEDFILE(ext) ;
    if (ffile->index == NULL){
      Rf_error("Failed to load bed index");
    }
    conf->bed = ffile->index;
  }

  if(cregion){
    conf->reg = cregion;
  }

  conf->output_fname = coutfn;

  // read the header of each file in the list and initialize data
  for (i = 0; i < n; ++i) {
    bam_hdr_t *h_tmp;
    data[i] = calloc(1, sizeof(mplp_aux_t));
    data[i]->fp = sam_open(cbampaths[i], "rb");

    if ( !data[i]->fp )
    {
      REprintf("failed to open %s: %s\n", cbampaths[i], strerror(errno));
      exit(EXIT_FAILURE);
    }

    if (conf->fai_fname) {
      if (hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
        REprintf("failed to process %s: %s\n",
                 cfapath, strerror(errno));
        exit(EXIT_FAILURE);
      }
    }

    data[i]->conf = conf;
    data[i]->ref = &mp_ref;

    h_tmp = sam_hdr_read(data[i]->fp);
    if ( !h_tmp ) {
      REprintf("fail to read the header of %s\n", cbampaths[i]);
      exit(EXIT_FAILURE);
    }

    if(conf->reg){
      hts_idx_t *idx = NULL;
      idx = sam_index_load(data[i]->fp, cbampaths[i]) ;

      if (idx == NULL) {
        REprintf("fail to load bamfile index for %s\n", cbampaths[i]);
        return 1;
      }

      if ( (data[i]->iter=sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
        REprintf("Fail to parse region '%s' with %s\n", conf->reg, cbampaths[i]);
        return 1;
      }

      if(i == 0){
        beg0 = data[i]->iter->beg;
        end0 = data[i]->iter->end;
        tid0 = data[i]->iter->tid;
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

  // return type void in rhtslib (htslib v 1.7), current v1.14 htslib returns int
  // enable overlap detection
  bam_mplp_init_overlaps(iter) ;

  // set max depth
  bam_mplp_set_maxcnt(iter, max_depth);

  if (!iter) {
    REprintf("issue with iterator");
    return 1;
  }

  pileup_fp = conf->output_fname? fopen(conf->output_fname, "w") : stdout;
  if (pileup_fp == NULL) {
    REprintf("Failed to write to %s: %s\n", conf->output_fname, strerror(errno));
    return 1;
  }

  pcounts *pc;
  pc = calloc(n, sizeof(pcounts));

  int min_1, min_2;
  min_1 = min_2 = 0;
  if(n == 1){
    min_1 = min_reads[0];
    min_2 = min_1;
  } else if (n == 2){
    min_1 = min_reads[0];
    min_2 = min_reads[1];
  }

  int n_iter = 0;
  int l;
  while ((l = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
      // check user interrupt
      // using a 2^k value (e.g. 256) can be 2-3x faster than say 1e6
      if (n_iter % 262144 == 0) R_CheckUserInterrupt();
      if (cregion && (pos < beg0 || pos >= end0)) continue; // out of the region requested
      mplp_get_ref(data[0], tid, &ref, &ref_len);

      if (tid < 0) break;
      if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, data[0]->h->target_name[tid], pos, pos+1)) continue;

      int pref_b, mref_b;
      pref_b = (ref && pos < ref_len)? ref[pos] : 'N' ;
      mref_b = comp_base[(unsigned char) pref_b];

      for(i = 0; i < n; ++i){
        memset(&pc[i], 0, sizeof(pcounts));
      }

      for (i = 0; i < n; ++i) {
        int j;

        if (n_plp[i] == 0) continue;

  	      // want to count plus and minus reads separately
  	      // iterate through reads that overlap position
        for (j = 0; j < n_plp[i]; ++j) {

          const bam_pileup1_t *p = plp[i] + j;

          // check base quality
          int bq = p->qpos < p->b->core.l_qseq
                   ? bam_get_qual(p->b)[p->qpos]
                   : 0;
          if (bq < min_baseQ) continue ;

          // skip indel and ref skip ;
          if(p->is_del || p->is_refskip) continue ;

          // consider adding a counter to track each error
          // to be used to exclude sites
          int check_artifacts = 0;
          if(check_artifacts){
            // check for splice in alignment nearby
            if(dist_to_splice(p->b, p->qpos, 4) != 0) continue;

            // check if indel event nearby
            if(dist_to_indel(p->b, p->qpos, 4) != 0) continue;

            // check if pos is within x dist from 5' end of read
            // qpos is 0-based
            if(trim_pos(p->b, p->qpos, 6)) continue;

          }

          // get read base
          int c = p->qpos < p->b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
          : 'N';

          int is_neg = 0;
          is_neg = bam_is_rev(p->b) ;

            // adjust based on library type and r1/r2 status
          int invert = 0;

          if(libtype == 1){

            if(!(p->b->core.flag&BAM_FPAIRED)){
              if(!(is_neg)) {
                invert = 1;
              }
            } else if (p->b->core.flag & (BAM_FREAD1)) {

              if(!(is_neg)) {
                invert = 1;
              }

            } else if (p->b->core.flag & (BAM_FREAD2)) {

              if(is_neg){
                invert = 1;
              }
            }

          } else if (libtype == 2){

            if(!(p->b->core.flag&BAM_FPAIRED)){
              if(is_neg) {
                invert = 1;
              }
            } else if (p->b->core.flag & (BAM_FREAD1)) {

              if(is_neg){
                invert = 1;
              }
            } else if (p->b->core.flag & (BAM_FREAD2)) {

              if(!(is_neg)) {
                invert = 1;
              }
            }
          }

  	      // NEED TO DEAL WITH UNSTRANDED
  	      // leave as positive strand for now?
          // } else {
          //   if(is_neg) {
  	      //     invert = 1;
          //   }
          // }

  	      // THIS COULD BE CLEANED UP SO LESS REPETITIVE
  	      // count reads that align to minus strand
          if(invert){
            c = comp_base[(unsigned char)c];

            pc[i].mtotal += 1;

            if(mref_b == c){
              pc[i].mnr += 1;
            } else {
              pc[i].mnv += 1;
            }

            switch(c) {
              case 'A':
                pc[i].mna += 1;
                break;
              case 'C':
                pc[i].mnc += 1;
                break;
              case 'G':
                pc[i].mng += 1;
                break;
              case 'T':
                pc[i].mnt += 1;
                break;
              default:
                pc[i].mnn += 1;
                break;
            }

          // count reads that align to plus strand
          } else {

            pc[i].ptotal += 1;

            if(pref_b == c){
              pc[i].pnr += 1;
            } else {
              pc[i].pnv += 1;
            }

            switch(c) {
              case 'A':
                pc[i].pna += 1;
                break;
              case 'C':
                pc[i].pnc += 1;
                break;
              case 'G':
                pc[i].png += 1;
                break;
              case 'T':
                pc[i].pnt += 1;
                break;
              default:
                pc[i].pnn += 1;
                break;
            }
          }
        }
      }
      // print out lines if pass depth criteria
      print_counts(pileup_fp, pc, n, min_1, min_2, h->target_name[tid], pos, pref_b, mref_b);
  }

  //clean up memory
  bam_mplp_destroy(iter);
  bam_hdr_destroy(h);

  for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
  free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);

  for (i = 0; i < n; ++i) {
    sam_close(data[i]->fp);
    if (data[i]->iter) hts_itr_destroy(data[i]->iter);
    free(data[i]);
  }
  free(data); free(plp); free(n_plp);
  free(pc);
  free(mp_ref.ref[0]);
  free(mp_ref.ref[1]);

  if (conf->fai) fai_destroy(conf->fai);

  // don't destroy index if passed from R
  if (conf->bed && Rf_isNull(ext)) bed_destroy(conf->bed);

  if (pileup_fp && conf->output_fname) fclose(pileup_fp);
  return 0;
}
