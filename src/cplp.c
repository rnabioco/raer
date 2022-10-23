#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/hts_log.h>
#include <bedidx.h>
#include "bedfile.h"
#include "utils.h"
#include "plp_data.h"

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

#include <time.h>

#include <Rinternals.h>

#define STRICT_R_HEADERS
#define R_NO_REMAP

//From https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

KHASH_SET_INIT_STR(rname)
typedef khash_t(rname) *rnhash_t;

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
  int min_mq, flag, min_baseQ, max_depth, all, multi, output_reads;
  read_qual_t read_qual;
  uint32_t keep_flag[2];
  char *reg, *fai_fname, *output_fname;
  faidx_t *fai;
  void *bed;
  rnhash_t rnames;
  strhash_t brhash;
  FILE *reads_fp;
  int argc;
  char **argv;
  char sep, empty;
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


/* from bam_plcmd.c */
static int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, hts_pos_t *ref_len) {

  mplp_ref_t *r = ma->ref;

  //REprintf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

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

  r->ref[0] = faidx_fetch_seq64(ma->conf->fai,
                              sam_hdr_tid2name(ma->h, r->ref_id[0]),
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

// read processing function for pileup
static int readaln(void *data, bam1_t *b) {
  char *ref;
  mplp_aux_t *g = (mplp_aux_t *)data;
  int ret, skip = 0;
  uint32_t test_flag;
  hts_pos_t ref_len;
  do {
    int has_ref;
    ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->h, b);
    if (ret < 0) break;

    // exclude unmapped reads
    if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) {
      skip = 1;
      continue;
    }

    // exclude reads in the "bad read hash"
    // should only occur if 1 bamfile is passed
    if (g->conf->brhash) {
      char* key;
      char c;
      size_t len = strlen(bam_get_qname(b));

      if(!((b)->core.flag&BAM_FPAIRED)){
        key = malloc(len + 1);
        if(!key) {
          Rf_error( "malloc failed\n");
        }
        strcpy(key, bam_get_qname(b));
        key[len] = '\0';
        c = '0';
      } else {
        key = malloc(len + 2 + 1);
        if(!key) {
          Rf_error( "malloc failed\n");
        }
        if((b)->core.flag&BAM_FREAD1){
          c = '1' ;
        } else {
          c = '2';
        }
        strcpy(key, bam_get_qname(b));
        key[len] = '_';
        key[len + 1] = c;
        key[len + 2] = '\0';
      }

      if (!key || kh_get(str, g->conf->brhash, key) != kh_end(g->conf->brhash)) {
        free(key);
        skip = 1;
        continue;
      }
      if(key) free(key);
    }

    // test required and filter flags
    test_flag = (g->conf->keep_flag[0] & ~b->core.flag) |
      (g->conf->keep_flag[1] & b->core.flag);
    if (~test_flag & 2047u){ skip = 1; continue;}

    // test overlap
    if (g->conf->bed && !g->conf->multi) {
      skip = !bed_overlap(g->conf->bed, sam_hdr_tid2name(g->h, b->core.tid), b->core.pos, bam_endpos(b));
      if (skip) continue;
    }

    if (g->conf->fai && b->core.tid >= 0) {
      has_ref = mplp_get_ref(g, b->core.tid, &ref, &ref_len);
      if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
        skip = 1;
        continue;
      }
    } else {
      has_ref = 0;
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

static void clear_rname_set(rnhash_t rnames)
{
  khint_t k;
  if(rnames){
    for (k = kh_begin(rnames); k < kh_end(rnames); ++k){
      if (kh_exist(rnames, k)){
        free((char*)kh_key(rnames, k));
      }
      kh_clear(rname, rnames);
    }
  }
}

// structs for holding counts
typedef struct {
  int total, nr, nv, na, nt, ng, nc, nn, nx;
  int ref_b;
  varhash_t var;
} counts;

typedef struct  {
  int pos;
  counts *pc;
  counts *mc;
} pcounts;

static void clear_counts(counts *p){
  p->na = p->nc = p->ng = p->nn = p->nr = p->nt = p->nv = p->total = p->nx = 0;
  p->ref_b = 0;
  clear_varhash_set(p->var);
}

static void clear_pcounts(pcounts *p){
  clear_counts(p->mc);
  clear_counts(p->pc);
  p->pos = 0;
}

// struct for event filter params
typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm, min_overhang;
} efilter;

static void get_var_string(varhash_t *vhash, char *out){
  const char *reg;
  int z = 1;
  khint_t k;

  if(kh_size(*vhash) > 0){
    for (k = kh_begin(*vhash); k < kh_end(*vhash); k++) {
      if (kh_exist(*vhash,k)) {
        reg = kh_key(*vhash,k);

        if(z == 1){
          strcpy(out, reg);
        } else {
          strcat(out, ",");
          strcat(out, reg);
        }
        z += 1;
      }
    }
  } else {
    strcpy(out, "-");
  }
}

static void print_counts(FILE *fp, counts *p, const char* ctig, int pos, int ref, int strand){
  char vout[12];
  get_var_string(&(p->var), vout);
  fprintf(fp,
          "%s\t%i\t%c\t%c\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
          ctig,
          pos + 1,
          strand,
          ref,
          vout,
          p->nr,
          p->nv,
          p->na,
          p->nt,
          p->nc,
          p->ng,
          p->nn,
          p->nx);
}

static void add_counts(PLP_DATA pd, int fi, counts *p, const char* ctig, int gpos, int ref, int strand){

  if(fi >= pd->nfiles || fi < 0){
    Rf_error("[raer internal] issue with file index %d", fi);
  }

  SEXP r = get_or_grow_PLP_DATA(pd, -1), s;
  int li;
  int idx = pd->icnt;
  char *buf;
  char vout[12];
  get_var_string(&(p->var), vout);

  for (li = 0; li < LENGTH(r); ++li) {
    if (R_NilValue == (s = VECTOR_ELT(r, li)))
      continue;
    switch (li) {
    case SEQNAME_IDX:
      buf = R_Calloc(strlen(ctig) + 1, char);
      if (!buf)
        Rf_error("add_plp_data: failed to allocate memory");
      strcpy(buf, ctig);
      pd->pdat[fi].seqnames[idx] = buf;
      break;
    case POS_IDX:
      pd->pdat[fi].pos[idx] = gpos + 1;
      break;
    case STRAND_IDX:
      buf = R_Calloc(2, char);
      buf[0] = (char)strand;
      pd->pdat[fi].strand[idx] = buf;
      break;
    case REF_IDX:
      buf = R_Calloc(2, char);
      buf[0] = (char)ref;
      pd->pdat[fi].ref[idx] = buf;
      break;
    case VAR_IDX:
      buf = R_Calloc(strlen(vout) + 1, char);
      strcpy(buf, vout);
      pd->pdat[fi].var[idx] = buf;
      break;
    case NREF_IDX:
      pd->pdat[fi].nref[idx] = p->nr;
      break;
    case NVAR_IDX:
      pd->pdat[fi].nvar[idx] = p->nv;
      break;
    case NA_IDX:
      pd->pdat[fi].na[idx] = p->na;
      break;
    case NT_IDX:
      pd->pdat[fi].nt[idx] = p->nt;
      break;
    case NC_IDX:
      pd->pdat[fi].nc[idx] = p->nc;
      break;
    case NG_IDX:
      pd->pdat[fi].ng[idx] = p->ng;
      break;
    case NN_IDX:
      pd->pdat[fi].nn[idx] = p->nn;
      break;
    case NX_IDX:
      pd->pdat[fi].nx[idx] = p->nx;
      break;

    default:
      Rf_error("[raer internal]: unhandled add_counts");
    break;
    }
  }
}


static void store_counts(PLP_DATA pd, pcounts *pc, int *only_variants, int n, int min_depth,
                         const char* ctig, const int gpos, int pref, int mref, int in_memory){
  int i;
  int pv = 0;
  int mv = 0;
  int write_p = 0;
  int write_m = 0;
  int write_only_v = 0;

  for(i = 0; i < n; ++i){
    if(only_variants[i]){
      write_only_v = 1;
    }
  }
  // write out if any samples are > min_depth and any samples have a variant
  // predicated on the true or false values in only_variants
  for(i = 0; i < n; ++i){
    if((pc + i)->pc->total >= min_depth){
      write_p = 1;
    }
    if((pc + i)->mc->total >= min_depth){
      write_m = 1;
    }
    if(write_only_v && only_variants[i]){
      if(kh_size((pc + i)->pc->var) > 0){
        pv = 1;
      }
      if(kh_size((pc + i)->mc->var) > 0){
        mv = 1;
      }
    }
  }
  if(write_only_v){
    write_p = pv && write_p;
    write_m = mv && write_m;
  }

  if(write_p){
    for(i = 0; i < n; ++i){
      if(in_memory){
        add_counts(pd, i, (pc + i)->pc, ctig, gpos, pref, '+');
      } else {
        print_counts(pd->fps[i], (pc + i)->pc, ctig, gpos, pref, '+') ;
      }
    }
    pd->icnt += 1;
  }
  if(write_m){
    for(i = 0; i < n; ++i){
      if(in_memory){
        add_counts(pd, i, (pc + i)->mc, ctig, gpos, mref, '-');
      } else {
        print_counts(pd->fps[i], (pc + i)->mc, ctig, gpos, mref, '-') ;
      }
    }
    pd->icnt += 1;
  }
}


static int write_fasta(bam1_t *b, mplp_conf_t *conf, const char* ref, char read, int pos){
  if(!conf->reads_fp){
    return -1;
  }
  char* seq;
  seq = get_read(b);
  fprintf(conf->reads_fp,
          ">%s_%c_%s_%i\n",
          bam_get_qname(b),
          read,
          ref,
          pos) ;
  fprintf(conf->reads_fp, "%s\n", seq);
  free(seq);
  return 0;
}

// add/check for readname (with 1 or 2 appended if paired) to hash table
// if unique write to bam file
static int write_reads(bam1_t *b,  mplp_conf_t *conf, const char* ref, const int pos){

  char* key;
  char c;
  size_t len = strlen(bam_get_qname(b));

  if(!((b)->core.flag&BAM_FPAIRED)){
    key = malloc(len + 1);
    if(!key) {
      Rf_error( "malloc failed\n");
      exit(EXIT_FAILURE);
    }
    strcpy(key, bam_get_qname(b));
    key[len] = '\0';
    c = '0';
  } else {
    key = malloc(len + 2 + 1);
    if(!key) {
      Rf_error( "malloc failed\n");
      exit(EXIT_FAILURE);
    }
    if((b)->core.flag&BAM_FREAD1){
      c = '1' ;
    } else {
      c = '2';
    }
    strcpy(key, bam_get_qname(b));
    key[len] = '_';
    key[len + 1] = c;
    key[len + 2] = '\0';
  }

  // try to add read to cache, check if already present
  // key will be freed upon cleanup of hashtable, i believe
  int rret;
  kh_put(rname, conf->rnames, key, &rret);
  if (rret == 1) {
    int ret = 0;
    ret = write_fasta(b, conf, ref, c, pos + 1);
    if(ret < 0) {
      Rf_error("writing read failed\n");
      exit(EXIT_FAILURE);
    }
  } else {
    free(key);
  }
  return 0;
}

static int count_one_record(bam1_t *b, pcounts *pc, mplp_conf_t *conf,
                            const char* ctig, int pos, int pref_b, int mref_b,
                            int read_base, int invert, int i){
  int hret = 0;

  char mvar[3];
  char pvar[3];
  mvar[2] = '\0';
  pvar[2] = '\0';

  if(invert){
    pc->mc->total += 1;

    // check if read base == ref base
    if(mref_b == read_base){
      pc->mc->nr += 1;
    } else {
      mvar[0] = mref_b;
      mvar[1] = read_base;

      // store variants as (AG, AT, AC, etc.) hash set
      char *var = strdup(mvar);
      kh_put(varhash, pc->mc->var, var, &hret);
      if (hret == 0) free(var);

      if(i == 0 && conf->output_reads){
        int wret = write_reads(b, conf, ctig, pos);
        if(wret != 0){
          REprintf( "writing mismatched reads failed\n");
          return -1;
        }
      }
      pc->mc->nv += 1;
    }
    switch(read_base) {
    case 'A':
      pc->mc->na += 1;
      break;
    case 'C':
      pc->mc->nc += 1;
      break;
    case 'G':
      pc->mc->ng += 1;
      break;
    case 'T':
      pc->mc->nt += 1;
      break;
    default:
      pc->mc->nn += 1;
    break;
    }

    // count reads from plus strand
  } else {
    pc->pc->total += 1;

    if(pref_b == read_base){
      pc->pc->nr += 1;
    } else {

      pvar[0] = pref_b;
      pvar[1] = read_base;

      char *var = strdup(pvar);
      kh_put(varhash, pc->pc->var, var, &hret);
      if (hret == 0) free(var);
      if(i == 0 && conf->output_reads){
        int wret = write_reads(b, conf, ctig, pos);
        if(wret != 0){
          REprintf( "writing mismatched reads failed\n");
          return -1;
        }
      }
      pc->pc->nv += 1;
    }
    switch(read_base) {
    case 'A':
      pc->pc->na += 1;
      break;
    case 'C':
      pc->pc->nc += 1;
      break;
    case 'G':
      pc->pc->ng += 1;
      break;
    case 'T':
      pc->pc->nt += 1;
      break;
    default:
      pc->pc->nn += 1;
    break;
    }
  }
  return 0;
}

/* return codes,
 * 0 = read passes, 
 * 1 = read fails filter should be counted as bad read
 * 2 = read fails do no count as bad read
 */
static int check_read_filters(const bam_pileup1_t *p, efilter *ef, int baq, int maq){

 // skip indel and ref skip ;
  if(p->is_del || p->is_refskip) return(2) ;

  // filter based on mapq as not able to filter per file in pileup
  if (p->b->core.qual < maq) return(1) ;

  // check base quality
  int bq = p->qpos < p->b->core.l_qseq
    ? bam_get_qual(p->b)[p->qpos]
    : 0;
  if (bq < baq) {
    if(bq == 0) { // overlapping read pairs get base qualities set to 0 by mplp 
      return(2);
    } else {
      return(1);
    } 
  } 
    
  // check if pos is within x dist from 5' end of read, qpos is 0-based
  if(trim_pos(p->b, p->qpos, ef->trim_5p_dist, ef->trim_3p_dist)) return(1);

  // check for splice in alignment nearby
  if(ef->splice_dist && dist_to_splice(p->b, p->qpos, ef->splice_dist) >= 0) return(1);

  // check if site in splice overhang and > min_overhang
  if(ef->min_overhang && check_splice_overhang(p->b, p->qpos, ef->min_overhang) > 0) return(1);

  // check if indel event nearby
  if(ef->indel_dist && dist_to_indel(p->b, p->qpos, ef->indel_dist) >= 0) return(1);

  return 0;
}

SEXP run_cpileup(const char** cbampaths,
                int n,
                char* cfapath,
                char* cregion,
                int in_mem,
                int multi_region_itr,
                char** coutfns,
                char* cbedfn,
                int min_reads,
                int max_depth,
                int min_baseQ,
                int* min_mapQ,
                int* libtype,
                int* b_flags,
                int* event_filters,
                int* only_keep_variants,
                char* reads_fn,
                char* mismatches,
                double* read_bqual_filter,
                SEXP ext) {

  hts_set_log_level(HTS_LOG_ERROR);

  mplp_aux_t **data;
  int i, tid, *n_plp, tid0 = 0, ret = 0;
  hts_pos_t pos, beg0 = 0, end0 = INT32_MAX, ref_len;

  if(min_baseQ < 0) min_baseQ = 0;
  if(!max_depth) max_depth = 10000;

  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  bam_mplp_t iter;
  bam_hdr_t *h = NULL; /* header of first file in input list */
  char *ref;

  SEXP result = PROTECT(pileup_result_init(n));
  PLP_DATA pd = init_PLP_DATA(result,  n);

  if(!in_mem) pd->fps = R_Calloc(n, FILE*);

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

  // load and index bed intervals or use pointer to index
  // optionally build a multi-region iterator
  if(cbedfn){
    conf->bed = bed_read(cbedfn);
    conf->multi = multi_region_itr;
  } else if (!Rf_isNull(ext)){
    _BED_FILE *ffile = BEDFILE(ext) ;
    if (ffile->index == NULL){
      Rf_error("Failed to load bed index");
    }
    conf->bed = ffile->index;
    conf->multi = multi_region_itr;
  }
  // single region for pileup
  if(cregion) conf->reg = cregion;

  if(mismatches){
    int chk = 0;
    if(n > 1){
      Rf_error("unable to exclude bad reads with multiple input files");
    }
    conf->brhash = kh_init(str);
    if (conf->brhash == NULL) {
      Rf_error("issue building hash");
    }
    chk = populate_lookup_from_file(conf->brhash, mismatches);
    if(chk < 0){
      Rf_error("issue building hash");
    }
  }

  // if multiple bam files, use minimum mapQ for initial filtering,
  // then later filter in pileup loop
  if(min_mapQ[0] >= 0) {
    int mq;
    mq = min_mapQ[0];
    for(int i = 0; i < n; ++i){
      if(min_mapQ[i] < mq){
        mq = min_mapQ[i];
      }
    }
    conf->min_mq = mq;
  }

  if(b_flags){
    conf->keep_flag[0] = b_flags[0];
    conf->keep_flag[1] = b_flags[1];
  }

  efilter *ef;
  ef = R_Calloc(1, efilter);
  ef->trim_5p_dist = event_filters[0];
  ef->trim_3p_dist = event_filters[1];
  ef->splice_dist = event_filters[2];
  ef->indel_dist = event_filters[3];
  ef->nmer = event_filters[4];
  ef->n_mm_type = event_filters[5];
  ef->n_mm = event_filters[6];
  ef->min_overhang = event_filters[7];

  if(read_bqual_filter){
    conf->read_qual.pct = read_bqual_filter[0];
    conf->read_qual.minq = (int)read_bqual_filter[1];
  }

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
        tid0 = data[i]->iter->tid;
      }
      hts_idx_destroy(idx);

    } else if (conf->bed && conf->multi){
      data[i]->idx = sam_index_load(data[i]->fp, cbampaths[i]) ;

      if (data[i]->idx == NULL) {
        Rf_error("fail to load bamfile index for %s\n", cbampaths[i]);
      }

      bed_unify(conf->bed);

      if (!conf->bed) { // index is unavailable or no regions have been specified
        Rf_error("No regions or BED file have been provided. Aborting.");
      }

      int regcount = 0;
      hts_reglist_t *reglist = bed_reglist(conf->bed, 0, &regcount);
      if (!reglist) {
        Rf_error("Region list is empty or could not be created. ");
      }

      data[i]->iter = sam_itr_regions(data[i]->idx, h_tmp, reglist, regcount);
      if(!data[i]->iter) {
        Rf_error("Multi-region iterator could not be created. Aborting.");
      }

      if(i == 0){
        beg0 = data[i]->iter->beg;
        end0 = data[i]->iter->end;
        tid0 = data[i]->iter->tid;
      }
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

  if (reads_fn){
    conf->output_reads = 1;
    conf->reads_fp = fopen(R_ExpandFileName(reads_fn), "w");
    if (!conf->reads_fp) {
      Rf_error("failed to open %s: %s\n",R_ExpandFileName(reads_fn), strerror(errno));
    }
    conf->rnames = kh_init(rname);
  }

  iter = bam_mplp_init(n, readaln, (void**)data);

  // return type void in rhtslib (htslib v 1.7), current v1.14 htslib returns int
  // enable overlap detection
  bam_mplp_init_overlaps(iter) ;

  // set max depth
  bam_mplp_set_maxcnt(iter, max_depth);

  if (!iter) {
    REprintf("issue with iterator");
    ret = -1;
    goto fail;
  }

  pcounts *plpc;
  plpc = R_Calloc(n, pcounts);
  for (i = 0; i < n; ++i) {
     if(!in_mem){
       // initialize output files
       pd->fps[i] = fopen(R_ExpandFileName(coutfns[i]), "w");
       if (pd->fps[i] == NULL) {
         REprintf("Failed to open file outputfile %s\n", R_ExpandFileName(coutfns[i]));
         ret = -1;
         goto fail;
       }
     }
     plpc[i].pc = R_Calloc(1, counts);
     plpc[i].mc = R_Calloc(1, counts);

     //initialize variant string set (AG, AT, TC) etc.
     if (plpc[i].pc->var == NULL) {
       plpc[i].pc->var = kh_init(varhash);
     }
     if (plpc[i].mc->var == NULL) {
       plpc[i].mc->var = kh_init(varhash);
     }
  }

  int last_tid = -1;
  int n_iter = 0;
  while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    if (cregion && (pos < beg0 || pos >= end0)) continue; // not in of single region requested

    mplp_get_ref(data[0], tid, &ref, &ref_len); // not in of single region requested
    if (tid < 0) break;
    
    // check user interrupt, using a 2^k value is 2-3x faster than say 1e6
    if (n_iter % 262144 == 0) {
      if(checkInterrupt()){
        goto fail;
      }
    }
    // ensure position is in requested intervals
    if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, sam_hdr_tid2name(h, tid), pos, pos+1)) continue;

    // if writing out reads with mismatches, clear out hash table of read names at each chromosome,
    if(tid != last_tid && conf->output_reads){
      if (kh_size(conf->rnames)) {
        clear_rname_set(conf->rnames);
      }
      last_tid = tid;
    }
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
                 in_mem);
  }

fail:
  if(in_mem) finish_PLP_DATA(pd);

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
    for (k = 0; k < kh_end(conf->brhash); ++k)
      if (kh_exist(conf->brhash, k)) free((char*)kh_key(conf->brhash, k));
      kh_destroy(str, conf->brhash);
  }

  UNPROTECT(1);

  if(ret < 0) Rf_error("error detected during pileup");

  if(in_mem){
    return pd->result;
  } else {
    return Rf_ScalarInteger(ret);
  }
}
