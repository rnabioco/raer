#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include "samtools/bedidx.h"
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

KHASH_SET_INIT_STR(varhash);
typedef khash_t(varhash) *varhash_t;

KHASH_SET_INIT_STR(rname);
typedef khash_t(rname) *rnhash_t;

KHASH_SET_INIT_STR(str)
  typedef khash_t(str) *strhash_t;

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

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

/* from samtools bam_fastq.c

* Reverse a string in place.
* From http://stackoverflow.com/questions/8534274/is-the-strrev-function-not-available-in-linux.
* Author Sumit-naik: http://stackoverflow.com/users/4590926/sumit-naik
*/
  static char *reverse(char *str)
  {
    int i = strlen(str)-1,j=0;
    char ch;
    while (i>j) {
      ch = str[i];
      str[i]= str[j];
      str[j] = ch;
      i--;
      j++;
    }
    return str;
  }


/* from samtools bam_view.c */

static int populate_lookup_from_file(strhash_t lookup, char *fn)
{
  FILE *fp;
  char buf[1024];
  int ret = 0;
  fp = fopen(fn, "r");
  if (fp == NULL) {
    Rf_error("failed to open \"%s\" for reading", fn);
    return -1;
  }

  while (ret != -1 && !feof(fp) && fscanf(fp, "%1023s", buf) > 0) {
    char *d = strdup(buf);
    if (d != NULL) {
      kh_put(str, lookup, d, &ret);
      if (ret == 0) free(d); /* Duplicate */
    } else {
      ret = -1;
    }
  }
  if (ferror(fp)) ret = -1;
  if (ret == -1) {
    Rf_error("failed to read \"%s\"", fn);
  }
  fclose(fp);
  return (ret != -1) ? 0 : -1;
}

/* from samtools bam_fastq.c */
/* return the read, reverse complemented if necessary */
static char *get_read(const bam1_t *rec)
{
  int len = rec->core.l_qseq + 1;
  char *read = calloc(1, len);
  char *seq = (char *)bam_get_seq(rec);
  int n;

  if (!read) return NULL;

  for (n=0; n < rec->core.l_qseq; n++) {
    if (rec->core.flag & BAM_FREVERSE) read[n] = seq_nt16_str[seq_comp_table[bam_seqi(seq,n)]];
    else                               read[n] = seq_nt16_str[bam_seqi(seq,n)];
  }
  if (rec->core.flag & BAM_FREVERSE) reverse(read);
  return read;
}

typedef struct {
  int minq;
  double pct;
} read_qual_t;

typedef struct {
  int min_mq, flag, min_baseQ, max_depth, all;
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

/* from bam_plcmd.c */
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

// read processing function for pileup
static int readaln(void *data, bam1_t *b) {
  mplp_aux_t *g = (mplp_aux_t *)data;
  int ret;
  int skip = 0;
  uint32_t test_flag;
  while (1) {
    ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->h, b);
    if (ret < 0) break;

    // exclude unmapped reads
    if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) continue;

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

      if (!key || kh_get(str, g->conf->brhash, key) != kh_end(g->conf->brhash)) {
        free(key);
        continue;
      }
      if(key) free(key);
    }

    // test required and filter flags
    test_flag = (g->conf->keep_flag[0] & ~b->core.flag) |
      (g->conf->keep_flag[1] & b->core.flag);
    if (~test_flag & 2047u)
      continue;

    // test overlap
    if (g->conf->bed) {
      skip = !bed_overlap(g->conf->bed, g->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
      if (skip) continue;
    }

    // check mapping quality
    if (b->core.qual < g->conf->min_mq) continue;

    // check for proper pair, if paired end
    if((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) continue ;

    // check if read quality is > x in at least y% of read
    if(g->conf->read_qual.pct &&
       g->conf->read_qual.minq &&
       !read_base_quality(b, g->conf->read_qual.pct, g->conf->read_qual.minq)) continue;

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
  varhash_t pvar;
  varhash_t mvar;
} pcounts;

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

static void clear_varhash_set(varhash_t vhash)
{
  khint_t k;
  if(vhash){
    for (k = kh_begin(vhash); k < kh_end(vhash); ++k){
      if (kh_exist(vhash, k)){
        free((char*)kh_key(vhash, k));
      }
      kh_clear(varhash, vhash);
    }
  }
}

static void clear_pcounts(pcounts *p){
  p->mna = p->mnc = p->mng = p->mnn = p->mnr = p->mnt = p->mnv = p->mtotal = 0;
  p->pna = p->pnc = p->png = p->pnn = p->pnr = p->pnt = p->pnv = p->ptotal = 0;
  p->mref_b = p->pref_b = p->pos = 0;
  clear_varhash_set(p->mvar);
  clear_varhash_set(p->pvar);
}

// struct for event filter params
typedef struct  {
  int nmer, splice_dist, indel_dist, trim_5p_dist, trim_3p_dist;
  int n_mm_type, n_mm;
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

static void print_plus_counts(FILE *fp, pcounts *p, char* ctig, int pos, int pref){

  char pout[12];
  get_var_string(&(p->pvar), pout);
  fprintf(fp,
          "%s\t%i\t%c\t%c\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
          ctig,
          pos + 1,
          '+',
          pref,
          pout,
          p->pnr,
          p->pnv,
          p->pna,
          p->pnt,
          p->pnc,
          p->png,
          p->pnn);
}

static void print_minus_counts(FILE *fp, pcounts *p, char* ctig, int pos, int mref){

  char mout[12];
  get_var_string(&(p->mvar), mout);
  fprintf(fp,
          "%s\t%i\t%c\t%c\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
          ctig,
          pos + 1,
          '-',
          mref,
          mout,
          p->mnr,
          p->mnv,
          p->mna,
          p->mnt,
          p->mnc,
          p->mng,
          p->mnn);
}

static void print_counts(FILE **fps, pcounts *pc, int *only_variants, int n, int min_depth,
                         char* ctig, int pos, int pref, int mref){
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
    if((pc + i)->ptotal >= min_depth){
      write_p = 1;
    }
    if((pc + i)->mtotal >= min_depth){
      write_m = 1;
    }
    if(write_only_v && only_variants[i]){
      if(kh_size((pc + i)->pvar) > 0){
        pv = 1;
      }
      if(kh_size((pc + i)->mvar) > 0){
        mv = 1;
      }
    }
  }

  if(write_only_v){
    write_p = pv && write_p;
    write_m = mv && write_m;
  }

  for(i = 0; i < n; ++i){
    if(write_p){
      print_plus_counts(fps[i], (pc + i), ctig, pos, pref) ;
    }
    if(write_m){
      print_minus_counts(fps[i], (pc + i), ctig, pos, mref) ;
    }
  }
}

int write_fasta(bam1_t *b, mplp_conf_t *conf, char* ref, char read, int pos){
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
int write_reads(bam1_t *b,  mplp_conf_t *conf, char* ref, int pos){

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

// Parse MD tag and enumerate types and number of mismatches
// Determine if read passes n_mis and n_types thresholds
//
//
// returns: -1 if no MD tag
//           0 if read passes filter
//           >0 with # of unique mismatches if doesn't pass
//
// based on cigar and MD parsing code from htsbox (written by Heng Li)
// https://github.com/lh3/htsbox/blob/ffc1e8ad4f61291c676a323ed833ade3ee681f7c/samview.c#L117

int parse_mismatches(bam1_t* b, int pos, int n_types, int n_mis){
  int ret = 0;
  if(n_types < 1 && n_mis < 1){
    return ret;
  }
  const uint32_t *cigar = bam_get_cigar(b);
  int k, x, y, c;

  uint8_t *pMD = 0;
  char *p;
  int nm = 0;

  char vars[3];
  vars[2] = '\0';
  varhash_t vh;
  vh = kh_init(varhash);

  if ((pMD = bam_aux_get(b, "MD")) != 0) {
    for (k = x = y = 0; k < b->core.n_cigar; ++k) {
        int op = bam_cigar_op(cigar[k]);
        int len = bam_cigar_oplen(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
          x += len, y += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)
          y += len;
    }

    p = bam_aux2Z(pMD);

    y = 0;
    while (isdigit(*p)) {
      y += strtol(p, &p, 10);
      if (*p == 0) {
        break;
      } else if (*p == '^') { // deletion
        ++p;
        while (isalpha(*p)) ++p;
      } else {
        while (isalpha(*p)) {
          if (y >= x) {
            y = -1;
            break;
          }
          c = y < b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(b), y)]
          : 'N';
          nm += 1;

          vars[0] = *p;
          vars[1] = c;
          char *var = strdup(vars);
          int rval = 0;
          kh_put(varhash, vh, var, &rval);
          if(rval == -1){
            Rf_error("issue tabulating variants per read at, %s", bam_get_qname(b));
          } else if (rval == 0){
            free(var);
          }
          ++y, ++p;
        }
        if (y == -1) break;
      }
    }

    if (x != y) {
      REprintf("inconsistent MD for read '%s' (%d != %d); ignore MD\n", bam_get_qname(b), x, y);
      ret = -1;
    }

    if(kh_size(vh) >= n_types && nm >= n_mis){
      ret = kh_size(vh);
    }
    clear_varhash_set(vh);
    kh_destroy(varhash, vh);
  }
  return ret;

}


int run_cpileup(const char** cbampaths,
                int n,
                char* cfapath,
                char* cregion,
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

  mplp_aux_t **data;
  int i, tid, *n_plp, tid0 = 0, ret = 0;
  int pos, beg0 = 0, end0 = INT32_MAX, ref_len;

  if(min_baseQ < 0) min_baseQ = 0;
  if(!max_depth) max_depth = 10000;

  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;
  bam_mplp_t iter;
  bam_hdr_t *h = NULL; /* header of first file in input list */
  char *ref;

  FILE **pileup_fps;
  pileup_fps = malloc(n * sizeof(FILE *));

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
  ef = calloc(1, sizeof(efilter));
  ef->trim_5p_dist = event_filters[0];
  ef->trim_3p_dist = event_filters[1];
  ef->splice_dist = event_filters[2];
  ef->indel_dist = event_filters[3];
  ef->nmer = event_filters[4];
  ef->n_mm_type = event_filters[5];
  ef->n_mm = event_filters[6];

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

  int write_mismatched_reads = 0;
  if (reads_fn){
    write_mismatched_reads = 1;
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
    Rf_error("issue with iterator");
    return -1;
  }

  pcounts *pc;
  pc = calloc(n, sizeof(pcounts));
  for (i = 0; i < n; ++i) {
    // initialize output files
    pileup_fps[i] = fopen(R_ExpandFileName(coutfns[i]), "w");
    if ( pileup_fps[i] == NULL) {
      Rf_error("Failed to open file outputfile %s\n", R_ExpandFileName(coutfns[i]));
    }

    // initialize variant string hash (AG, AT, TC) etc.
    if (pc[i].mvar == NULL) {
      pc[i].mvar = kh_init(varhash);
    }
    if (pc[i].pvar == NULL) {
      pc[i].pvar = kh_init(varhash);
    }
  }

  char mvar[3];
  char pvar[3];
  mvar[2] = '\0';
  pvar[2] = '\0';

  int hret ;

  int last_tid = -1;
  int n_iter = 0;
  while ((ret = bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
      // check user interrupt
      // using a 2^k value (e.g. 256) can be 2-3x faster than say 1e6
      if (n_iter % 262144 == 0) {
        if(checkInterrupt()){
          goto fail;
        }
      }

      if (cregion && (pos < beg0 || pos >= end0)) continue; // out of the region requested
      mplp_get_ref(data[0], tid, &ref, &ref_len);

      if (tid < 0) break;
      if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, data[0]->h->target_name[tid], pos, pos+1)) continue;

      // clear out hash table of read names each chromosome
      if(tid != last_tid && write_mismatched_reads){
        if (kh_size(conf->rnames)) {
          clear_rname_set(conf->rnames);
        }
        last_tid = tid;
      }

      int pref_b, mref_b;
      pref_b = (ref && pos < ref_len)? ref[pos] : 'N' ;
      pref_b = toupper(pref_b);
      mref_b = comp_base[(unsigned char) pref_b];

      for(i = 0; i < n; ++i){
        clear_pcounts(&pc[i]);
      }

      // check if site is in a homopolymer
      if(ef->nmer > 0 && check_simple_repeat(&ref, &ref_len, pos, ef->nmer)) continue;

      int pass_reads = 0;
      for(i = 0; i < n; ++i){
        // fail early if read count less than min_reads
        if (n_plp[i] >= min_reads) {
          pass_reads = 1;
        }
      }
      if(!pass_reads){
        continue;
      }

      for (i = 0; i < n; ++i) {
        int j;

  	    // iterate through reads that overlap position
        for (j = 0; j < n_plp[i]; ++j) {

          const bam_pileup1_t *p = plp[i] + j;

          // filter based on mapq as not able to filter per file in pileup
          if (p->b->core.qual < min_mapQ[i]) continue ;

          // check base quality
          int bq = p->qpos < p->b->core.l_qseq
                   ? bam_get_qual(p->b)[p->qpos]
                   : 0;
          // note that overlap detection can fail if one mate has a cigar with
          // splicing events
          // this was fixed in htslib #802 and release v1.10
          if (bq < min_baseQ) continue ;

          // skip indel and ref skip ;
          if(p->is_del || p->is_refskip) continue ;

          // consider adding a counter to track each error
          // check if pos is within x dist from 5' end of read
          // qpos is 0-based
          if(trim_pos(p->b, p->qpos, ef->trim_5p_dist, ef->trim_3p_dist)) continue;

          // check for splice in alignment nearby
          if(ef->splice_dist && dist_to_splice(p->b, p->qpos, ef->splice_dist) >= 0) continue;

          // check if indel event nearby
          if(ef->indel_dist && dist_to_indel(p->b, p->qpos, ef->indel_dist) >= 0) continue;

          // get read base
          int c = p->qpos < p->b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
          : 'N';

          int is_neg = 0;
          is_neg = bam_is_rev(p->b) ;

            // adjust based on library type and r1/r2 status
          int invert = 0;

          if(libtype[i] == 1){

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

          } else if (libtype[i] == 2){

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
          } else if (libtype[i] == 3){
            if(is_neg){
              invert = 1;
            }
          }
          // if libtype == 0, do nothing, report as plus strand

  	      // count reads that align to minus strand
          if(invert){
            c = (char)comp_base[(unsigned char)c];
            mvar[0] = mref_b;
            pc[i].mtotal += 1;

            if(mref_b == c){
              pc[i].mnr += 1;
            } else {
              // store variants in khash
              mvar[1] = c;
              int m;
              // drop read if >= mismatch different types and at least n_mm mismatches
              m = parse_mismatches(p->b, pos, ef->n_mm_type, ef->n_mm);
              if(m > 0){
                // exclude read
                continue;
              }

              char *var = strdup(mvar);
              kh_put(varhash, pc[i].mvar, var, &hret);
              if (hret == 0) free(var);
              if(i == 0 && write_mismatched_reads){
                int wret = write_reads(p->b, conf, h->target_name[tid], pos);
                if(wret != 0){
                  REprintf( "writing mismatched reads failed\n");
                  ret = -1;
                  goto fail;
                }
              }

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
            pvar[0] = pref_b;
            pc[i].ptotal += 1;

            if(pref_b == c){
              pc[i].pnr += 1;
            } else {
              pvar[1] = c;
              int m;
              // drop read if >= mismatch different types and at least n_mm mismatches
              m = parse_mismatches(p->b, pos, ef->n_mm_type, ef->n_mm);
              if(m > 0){
                // exclude read
                continue;
              }

              char *var = strdup(pvar);
              kh_put(varhash, pc[i].pvar, var, &hret);
              if (hret == 0) free(var);
              if(i == 0 && write_mismatched_reads){
                int wret = write_reads(p->b, conf, h->target_name[tid], pos);
                if(wret != 0){
                  REprintf( "writing mismatched reads failed\n");
                  ret = -1;
                  goto fail;
                }
              }

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
      print_counts(pileup_fps, pc, only_keep_variants, n, min_reads, h->target_name[tid], pos, pref_b, mref_b);
  }

fail:
  //clean up memory
  bam_mplp_destroy(iter);
  bam_hdr_destroy(h);

  for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
  free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);

  for (i = 0; i < n; ++i) {
    sam_close(data[i]->fp);
    if (data[i]->iter) hts_itr_destroy(data[i]->iter);
    free(data[i]);
    clear_pcounts(&pc[i]);
    kh_destroy(varhash, pc[i].mvar);
    kh_destroy(varhash, pc[i].pvar);
    fclose(pileup_fps[i]);
  }
  free(data); free(plp); free(n_plp);
  free(pc); free(ef);
  free(mp_ref.ref[0]);
  free(mp_ref.ref[1]);
  free(pileup_fps);
  if (conf->fai) fai_destroy(conf->fai);

  // don't destroy index if passed from R
  if (conf->bed && Rf_isNull(ext)) bed_destroy(conf->bed);

  if (write_mismatched_reads) {
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

  return ret;
}
