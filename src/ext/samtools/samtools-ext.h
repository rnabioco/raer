#ifndef RAER_SAMTOOLS_EXT_H
#define RAER_SAMTOOLS_EXT_H

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include "../../plp_utils.h"

typedef struct {
  char *ref[2];
  int ref_id[2];
  hts_pos_t ref_len[2];
} mplp_ref_t;

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

int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, hts_pos_t *ref_len);

#endif
