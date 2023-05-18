#include <htslib/sam.h>
#include <htslib/khash.h>
#include <ctype.h>
#include <Rinternals.h>
#include "../../plp_utils.h"
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

int parse_mismatches(bam1_t* b, int n_types, int n_mis){
  int ret = 0;
  if(n_types < 1 && n_mis < 1){
    return ret;
  }
  const uint32_t *cigar = bam_get_cigar(b);
  int k, x, y, c;

  uint8_t *pMD = 0;
  char *p;
  int nm = 0;

  char vars[2];
  vars[1] = '\0';
  str2intmap_t vh;
  vh = kh_init(str2intmap);

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

          vars[0] = c;
          char *var = strdup(vars);
          int rval = 0;
          kh_put(str2intmap, vh, var, &rval);
          if(rval == -1){
            REprintf("[raer internal] issue tabulating variants per read at, %s\n",
                     bam_get_qname(b));
            clear_str2int_hashmap(vh);
            kh_destroy(str2intmap, vh);
            return -1;
          } else if (rval == 0){
            free(var);
          }
          ++y, ++p;
        }
        if (y == -1) break;
      }
    }

    if (x != y) {
      REprintf("[raer internal] inconsistent MD for read '%s' (%d != %d); ignore MD\n",
               bam_get_qname(b), x, y);
      ret = -1;
    } else {
      if(kh_size(vh) >= n_types && nm >= n_mis){
        ret = kh_size(vh);
      }
    }
  }
  clear_str2int_hashmap(vh);
  kh_destroy(str2intmap, vh);
  return ret;
}

