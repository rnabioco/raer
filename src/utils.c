#include <Rinternals.h>
#include "htslib/sam.h"

// https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
SEXP isnull(SEXP pointer) {
  return ScalarLogical(!R_ExternalPtrAddr(pointer));
}

int check_simple_repeat(char** ref, int* ref_len, int pos, int nmer){
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
  for (int i = 0;i < n_pos; i++) {
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

// check if query position is within dist from 5' end of read
int trim_pos(bam1_t* b, int pos, int dist){
  // pos is 0-based
  if(!(b->core.flag&BAM_FPAIRED) || b->core.flag&(BAM_FREAD1)){
    if(pos < dist){
      return 1;
    }
  } else if (b->core.flag&(BAM_FREAD2)) {
    if((b->core.l_qseq - pos) <= dist){
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
  for (i = 0; i < n_cigar; i++){
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

// return -1 if no indel found
int dist_to_indel(bam1_t* b, int pos, int dist){
  int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);

  int i, c, read_pos, indel_start, indel_end, sp, ep, ldist, rdist, idist;
  read_pos = indel_start = indel_end = 0;
  sp = pos - dist;
  ep = pos + dist;

  for (i = 0; i < n_cigar; i++){
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
