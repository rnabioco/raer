#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "bedidx.h"
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
    if(!aux_info && n_errors < 5) {
      REprintf("Warning: missing tag info in reads, e.g: %s", bam_get_qname(b)) ;
      n_errors += 1 ;
      continue ;
    }
    int64_t mmap_tag = bam_aux2i(aux_info);
    if(mmap_tag > 1) continue ;

    break;
  }

  return ret;
}


int run_cpileup(char* cbampath,
                char* cfapath,
                char* cregion,
                char* coutfn,
                char* cbedfn) {

  int n = 1;
  mplp_aux_t **data;
  int i, tid, *n_plp, tid0 = 0, max_depth;
  int pos, beg0 = 0, end0 = INT32_MAX, ref_len;
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

  // set this up in a loop
  data[0] = calloc(1, sizeof(mplp_aux_t));
  data[0]->fp = sam_open(cbampath, "rb");
  data[0]->h = sam_hdr_read(data[0]->fp);

  // read the header of each file in the list and initialize data
  if (cfapath) {
    if (hts_set_fai_filename(data[0]->fp, cfapath) != 0) {
      REprintf("failed to process %s: %s\n",
               cfapath, strerror(errno));
      exit(EXIT_FAILURE);
    }
  }
  if(cregion){
    hts_idx_t *idx = NULL;
    idx = sam_index_load(data[0]->fp, cbampath) ;

    if (idx == NULL) {
      REprintf("fail to load bamfile index for %s\n", cbampath);
      return 1;
    }

    if ( (data[0]->iter=sam_itr_querys(idx, data[0]->h, cregion)) == 0) {
      REprintf("Fail to parse region '%s' with %s\n", cregion, cbampath);
      return 1;
    }
    beg0 = data[0]->iter->beg;
    end0 = data[0]->iter->end;
    tid0 = data[0]->iter->tid;
    hts_idx_destroy(idx);

  } else {
    data[0]->iter = NULL;
  }

  if(cfapath){
    conf->fai_fname = cfapath;
    conf->fai = fai_load(conf->fai_fname);
  }

  if(cbedfn){
    conf->bed = bed_read(cbedfn) ;
  }

  conf->output_fname = coutfn;

  data[0]->conf = conf;
  data[0]->ref = &mp_ref;

  int status ;

  // init int array, if multiple files, set to n files
  int max_count = 10000;
  iter = bam_mplp_init(1, readaln, (void**)data);

  // return type void in rhtslib (htslib v 1.7), current v1.14 htslib returns int
  // enable overlap detection
  bam_mplp_init_overlaps(iter) ;

  // set max depth
  bam_mplp_set_maxcnt(iter, max_count);

  if (!iter) {
    REprintf("issue with iterator");
    return 1;
  }
  int min_baseQ = 20;
  int n_records = 0;
  int min_depth = 20;

  pileup_fp = conf->output_fname? fopen(conf->output_fname, "w") : stdout;
  if (pileup_fp == NULL) {
    REprintf("Failed to write to %s: %s\n", conf->output_fname, strerror(errno));
    return 1;
  }
  while ((n = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
      if (cregion && (pos < beg0 || pos >= end0)) continue; // out of the region requested
      mplp_get_ref(data[0], tid, &ref, &ref_len);
      if (tid < 0) break;
    //  if (n_records > 10) break ;
    //  Rprintf("n_reads_per_col: %i\n", n_plp[0]);

      if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, data[0]->h->target_name[tid], pos, pos+1)) continue;
      if (n_plp[0] > 1){

        int j;
        int ref_b;
        int total, nr, nv, na, nt, ng, nc, nn;
        total = nr = nv = na = nt = ng = nc = nn = 0;
        char strand = '+';

        for (j = 0; j < n_plp[0]; ++j) {
          const bam_pileup1_t *p = plp[0] + j;

          // check base quality
          int bq = p->qpos < p->b->core.l_qseq
            ? bam_get_qual(p->b)[p->qpos]
          : 0;
          if (bq < min_baseQ) continue ;

          // skip indel and ref skip ;
          if(p->is_del || p->is_refskip) continue ;

          //
          int c = p->qpos < p->b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
          : 'N';
          ref_b = (ref && pos < ref_len)? ref[pos] : 'N' ;

          // todo: adjust based on library type and r1/r2 status
          if(bam_is_rev(p->b)){
            c = comp_base[(unsigned char)c];
            ref_b = comp_base[(unsigned char)ref_b];
            strand = '-';
          }

          total += 1;
          if(ref_b == c){
            nr += 1;
          } else {
            nv += 1;
          }

          switch(c) {
            case 'A':
              na += 1;
              break;
            case 'C':
              nc += 1;
              break;
            case 'G':
              ng += 1;
              break;
            case 'T':
              nt += 1;
              break;
            default:
              nn += 1;
              break;
          }

        //  Rprintf("ref_nt: %c\n", ref_b);
        //  Rprintf("qname: %s\n", bam_get_qname(p->b));
        //  Rprintf("query_nt: %c\n", seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]);
        //  Rprintf("bq: %i\n", bq);
          n_records += 1 ;
        }
        if(total > min_depth){
          fprintf(pileup_fp,
                  "%s\t%i\t%c\t%c\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
                  data[0]->h->target_name[tid],
                  pos + 1,
                  strand,
                  ref_b,
                  nr,
                  nv,
                  na,
                  nt,
                  nc,
                  ng,
                  nn);
        }
      }
  }
  bam_mplp_destroy(iter);

  bam_hdr_destroy(data[0]->h);
  sam_close(data[0]->fp);
  if (pileup_fp && conf->output_fname) fclose(pileup_fp);
  return 0;
}
