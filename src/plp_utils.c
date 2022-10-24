#include <Rinternals.h>
#include <htslib/sam.h>
#include <htslib/khash.h>

#include <ctype.h>

#include "plp_utils.h"

// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer) {
  return ScalarLogical(!R_ExternalPtrAddr(pointer));
}

// Based on pysam code Cython code
// https://github.com/pysam-developers/pysam/blob/ef06e42ce98e4c81972a448ddab62289bf3ad22d/pysam/libcalignedsegment.pyx#L501

int query_start(bam1_t *b){
  const uint32_t *cigar = bam_get_cigar(b);
  int start_offset = 0;
  int n_cigar = b->core.n_cigar;
  int i, op;

  for (i = 0; i <  n_cigar; ++i){
    op = cigar[i] & BAM_CIGAR_MASK;
    if (op == BAM_CHARD_CLIP){
      if (start_offset != 0 && start_offset != b->core.l_qseq){
        Rf_error("Invalid clipping in CIGAR string");
      }
    } else if (op == BAM_CSOFT_CLIP){
      start_offset += cigar[i] >> BAM_CIGAR_SHIFT;
    } else {
      break;
    }
  }
  return start_offset;
}

int query_end(bam1_t *b){
  const uint32_t *cigar = bam_get_cigar(b);
  int n_cigar = b->core.n_cigar;
  unsigned int i;
  int op;
  int end_offset = b->core.l_qseq;

  if(end_offset == 0){
    Rf_error("SEQ record missing from BAM file");
  }

  for (i = n_cigar; i-- > 0;){
    op = cigar[i] & BAM_CIGAR_MASK;
    if (op == BAM_CHARD_CLIP){
      if (end_offset != 0 && end_offset != b->core.l_qseq){
        Rf_error("Invalid clipping in CIGAR string");
      }
    } else if (op == BAM_CSOFT_CLIP){
      end_offset += cigar[i] >> BAM_CIGAR_SHIFT;
    } else {
      break;
    }
  }

  return end_offset;
}

int check_simple_repeat(char** ref, hts_pos_t* ref_len, int pos, int nmer){
  int start, n_pos;
  n_pos = (nmer * 2) - 1;

  if(n_pos < 1 || nmer < 1 || *ref_len < nmer || pos < 0){
    return 0;
  }
  start = pos - nmer + 1;
  if(start < 0){
    n_pos = n_pos + start;
    start = 0;
  }
  n_pos = ((n_pos + start) > *ref_len) ? (*ref_len - start): n_pos;

  // check for homopolymer using run-length-encoding approach
  for (int i = 0;i < n_pos; ++i) {
    // Count occurrences of current character
    int count = 1;
    while (i < n_pos - 1 && (*ref)[start + i] == (*ref)[start + i + 1]) {
      count++;
      i++;
    }
    if(count >= nmer){
      return 1;
    }
  }
  return 0;
}

// check if query position is within dist from 5' or 3' end of alignment
// pos = query position
int trim_pos(bam1_t* b, int pos, int dist_5p, int dist_3p){
  // pos is 0-based query position
  // need to adjust to trim based on alignment start/end
  int qs, qe;
  qs = query_start(b);
  qe = query_end(b);

  if(!(b->core.flag&BAM_FREVERSE)){
    if(pos < (dist_5p + qs) || (qe - pos) <= dist_3p){
      return 1;
    }
  } else if (b->core.flag&(BAM_FREVERSE)) {
    if((qe - pos) <= dist_5p || pos < (dist_3p + qs)){
      return 1;
    }
  } else {
    Rf_error("don't believe should happen in trim_pos");
  }
  return 0;
}

// return -1 if no splice found
int dist_to_splice(bam1_t* b, int pos, int dist){
  int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);

  int i, c, ip, sp, ep, idist;
  ip = 0;
  sp = pos - dist;
  ep = pos + dist;
  for (i = 0; i < n_cigar; ++i){
    c = bam_cigar_op(cigar[i]);
    // is a M, I, S, =, or other query consuming operation
    if (bam_cigar_type(c)&1){
      ip += bam_cigar_oplen(cigar[i]);
      continue;
    }

    if(c == BAM_CREF_SKIP){
      if (ip >= sp && ip <= ep){
        idist = ip - pos;
        idist = (idist < 0) ? -idist : idist;
        return idist;
      }
    }
  }
  return -1;
}

// return
// -1 if overhand >= dist,
// 0 if no flanking splice
// or overlang length if less than dist
int check_splice_overhang(bam1_t* b, int pos, int dist){
  int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);

  int i, c, ip, cl, p_op, r;
  ip = 0;
  p_op = -1; // CIGAR macros are from 0-9,
  for (i = 0; i < n_cigar; ++i){
    c = bam_cigar_op(cigar[i]);
    cl = bam_cigar_oplen(cigar[i]);
    if(c == BAM_CMATCH && pos >= ip && pos <= (cl + ip)) {
      // site is in the 3' exon
      if(i > 0 && p_op == BAM_CREF_SKIP){
        r = cl >= dist ? -1 : cl;
        return r;
      }
      // site is in the 5' exon
      if(i < n_cigar){
        ++i;
        c = bam_cigar_op(cigar[i]);
        if(c == BAM_CREF_SKIP){
          r = cl >= dist ? -1 : cl;
          return r;
        }
      }
      // site not flanked by splice
      return 0;
    }
    p_op = c;
    if (bam_cigar_type(c)&1){
      ip += cl;
      continue;
    }
  }
  Rf_error("site not found in read: %s %i %i %i", bam_get_qname(b), pos, p_op, cl);
}

// return -1 if no indel found
int dist_to_indel(bam1_t* b, int pos, int dist){
  int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);

  int i, c, read_pos, indel_start, indel_end, sp, ep, ldist, rdist, idist;
  read_pos = indel_start = indel_end = 0;
  sp = pos - dist;
  ep = pos + dist;

  for (i = 0; i < n_cigar; ++i){
    c = bam_cigar_op(cigar[i]);

    // is a M, I, S, =, or other query consuming operation
    if (bam_cigar_type(c)&1){
      if(c == BAM_CINS){
        indel_start = read_pos;
        indel_end = read_pos + bam_cigar_oplen(cigar[i]);
        //check range overlap
        if (sp <= indel_end && indel_start <= ep){
          ldist = pos - indel_end;
          rdist = indel_start - pos;
          idist = (ldist > rdist) ? ldist : rdist;
          return idist;
        }
      }
      read_pos += bam_cigar_oplen(cigar[i]);
    }

    if(c == BAM_CDEL){
      if (read_pos >= sp && read_pos <= ep){
        idist = read_pos - pos;
        idist = (idist < 0) ? -idist : idist;
        return idist;
      }
    }
  }
  return -1;
}

/* return 1 if read passes quality thresholds
 * otherwise return 0 */
int read_base_quality(bam1_t* b, float pc, int mq){
  int c, i, ret, lq;
  ret = lq = 0;
  double qual_pct;

  if(!b->core.l_qseq){
    return 0;
  }

  for(i = 0; i < b->core.l_qseq; ++i){
    c = bam_get_qual(b)[i];
    if(c < mq){
      lq += 1;
    }
  }
  qual_pct = (double)lq / b->core.l_qseq;
  ret = qual_pct <= pc ? 1: 0;
  return ret;
}

// return 0 (plus) 1 (negative) -1 (error)
int invert_read_orientation(bam1_t* b, int libtype) {

  int is_neg = 0;
  is_neg = bam_is_rev(b) ;

  int invert = 0;
  if(libtype == 0){
    invert = 0;
  } else if (libtype == 1){
    if(!(b->core.flag&BAM_FPAIRED)){
      if(!(is_neg)) {
        invert = 1;
      }
    } else if (b->core.flag & (BAM_FREAD1)) {
      if(!(is_neg)) {
        invert = 1;
      }
    } else if (b->core.flag & (BAM_FREAD2)) {
      if(is_neg){
        invert = 1;
      }
    }
  } else if (libtype == 2){
    if(!(b->core.flag&BAM_FPAIRED)){
      if(is_neg) {
        invert = 1;
      }
    } else if (b->core.flag & (BAM_FREAD1)) {
      if(is_neg){
        invert = 1;
      }
    } else if (b->core.flag & (BAM_FREAD2)) {
      if(!(is_neg)) {
        invert = 1;
      }
    }
  } else if (libtype == 3){
    if(is_neg){
      invert = 1;
    }
  } else {
    invert = -1;
  }
  return invert;
}

void clear_varhash_set(varhash_t vhash)
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

int parse_mismatches(bam1_t* b,const int pos, int n_types, int n_mis){
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


/* from samtools bam_fastq.c

 * Reverse a string in place.
 * From http://stackoverflow.com/questions/8534274/is-the-strrev-function-not-available-in-linux.
 * Author Sumit-naik: http://stackoverflow.com/users/4590926/sumit-naik
 */
char *reverse(char *str)
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

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

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
int populate_lookup_from_file(strhash_t lookup, char *fn)
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
