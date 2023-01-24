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


char* get_aux_ztag(bam1_t *b, const char tag[2]){
  char* str;
  uint8_t* val = bam_aux_get(b, tag) ;
  if(val){
    str = bam_aux2Z(val);
  } else {
    str = NULL;
  }
  return(str);
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

unsigned char comp_base[256] = {
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
