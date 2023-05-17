#include <htslib/sam.h>
#include <htslib/faidx.h>
#include "../../plp_utils.h"
#include "samtools-ext.h"
#include <htslib/hts.h>

/* from bam_plcmd.c */
int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, hts_pos_t *ref_len) {

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

static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

/* from samtools bam_fastq.c */
/* return the read, reverse complemented if necessary */
char *get_read(const bam1_t *rec)
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

/* from samtools bam_view.c */
int populate_lookup_from_file(strset_t lookup, char *fn)
{
  FILE *fp;
  char buf[1024];
  int ret = 0;
  fp = fopen(fn, "r");
  if (fp == NULL) {
    REprintf("[raer internal] failed to open \"%s\" for reading", fn);
    return -1;
  }

  while (ret != -1 && !feof(fp) && fscanf(fp, "%1023s", buf) > 0) {
    char *d = strdup(buf);
    if (d != NULL) {
      kh_put(strset, lookup, d, &ret);
      if (ret == 0) free(d); /* Duplicate */
    } else {
      ret = -1;
    }
  }
  if (ferror(fp)) ret = -1;
  fclose(fp);
  if (ret == -1) {
    REprintf("[raer internal] failed to read \"%s\"", fn);
  }
  return (ret != -1) ? 0 : -1;
}
