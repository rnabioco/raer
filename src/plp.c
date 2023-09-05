#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/hts_log.h>

#include "regfile.h"
#include "plp_utils.h"
#include "plp_data.h"
#include "ext/bcftools/bcftools-ext.h"
#include "ext/samtools/samtools-ext.h"
#include "ext/htsbox/htsbox-ext.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define NBASE_POS 100
#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

// read processing function for pileup

/*! @function
 @abstract  Read processing function for pileup, selects reads to populate pileup

 @param  data     pointer to mplp_aux_t configuration
 @param  bam1_t   pointer to alignment struct

 @return
 0 if read passes filters
 1 otherwise
 */
static int readaln(void* data, bam1_t* b) {
    mplp_aux_t* g = (mplp_aux_t*)data;
    int ret, skip = 0;
    uint32_t test_flag;
    do {
        ret = g->iter? sam_itr_next(g->fp, g->iter, b) : sam_read1(g->fp, g->h, b);
        if (ret < 0) break;

        // exclude unmapped reads
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) {
            skip = 1;
            continue;
        }

        // test required and filter flags
        test_flag = (g->conf->keep_flag[0] & ~b->core.flag) |
            (g->conf->keep_flag[1] & b->core.flag);
        if (~test_flag & 2047u) {
            skip = 1;
            continue;
        }

        // test overlap
        if (g->conf->reg_idx) {
            skip = !regidx_overlap(g->conf->reg_idx,
                                   sam_hdr_tid2name(g->h, b->core.tid),
                                   b->core.pos,
                                   bam_endpos(b),
                                   g->conf->reg_itr);
            if (skip) continue;
        }

        skip = 0;
        // check mapping quality
        if (b->core.qual < g->conf->min_global_mq) {
            skip = 1;
        } else if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) {
            skip = 1;
            // check if read quality is > x in at least y% of read
        } else if (g->conf->read_qual.pct &&
            g->conf->read_qual.minq &&
            !read_base_quality(b, g->conf->read_qual.pct, g->conf->read_qual.minq)) {
            skip = 1;
        }
    } while (skip);
    return ret;
}

/*! @typedef
 @abstract Structure for storing count data for  bases
 @field var variants found at each site
 @field umi hashset with umi sequences at each site

 */
typedef struct {
    int total, nr, nv, na, nt, ng, nc, nn, nx;
    int ref_b;
    str2intmap_t var;
    strset_t umi;
} counts;

typedef struct  {
    int pos;
    counts* pc;
    counts* mc;
} pcounts;

/*! @typedef
 @abstract Structure for storing count data ignoring library-type for calculating
 strand bias
 */
typedef struct {
    int ref_fwd;
    int ref_rev;
    int alt_fwd;
    int alt_rev;
} strcounts_t;

/*! @typedef
 @abstract Structure for storing per-site data for all samples, used for assessing variant
 biases, such as strand bias
 */
typedef struct {
    int* p_ref_pos; // base position within read, (plus strand) scaled to 100 positions
    int* p_alt_pos; // variant base position within read, (plus strand) scaled to 100 positions
    int* m_ref_pos;
    int* m_alt_pos;
    int p_has_var; // 1 if any sample had a variant
    int m_has_var; // 1 if any sample had a variant
    strcounts_t* s; // counts for strandedness test
} pall_counts;

/*! @function
 @abstract  reset data to 0
 */
static void clear_pall_counts(pall_counts* p) {
    if (p->p_ref_pos) memset(p->p_ref_pos, 0, sizeof(int) * NBASE_POS);
    if (p->p_alt_pos) memset(p->p_alt_pos, 0, sizeof(int) * NBASE_POS);
    if (p->m_ref_pos) memset(p->m_ref_pos, 0, sizeof(int) * NBASE_POS);
    if (p->m_alt_pos) memset(p->m_alt_pos, 0, sizeof(int) * NBASE_POS);
    if (p->s) memset(p->s, 0, sizeof(strcounts_t));
    p->p_has_var = p->m_has_var = 0;
}

/*! @function
 @abstract  allocate data
 */
static pall_counts* init_pall_counts() {
    pall_counts* pall = R_Calloc(1, pall_counts);
    pall->p_ref_pos   = R_Calloc(NBASE_POS, int);
    pall->p_alt_pos   = R_Calloc(NBASE_POS, int);
    pall->m_ref_pos   = R_Calloc(NBASE_POS, int);
    pall->m_alt_pos   = R_Calloc(NBASE_POS, int);
    pall->s           = R_Calloc(1, strcounts_t);
    pall->p_has_var   = pall->m_has_var = 0;
    return (pall);
}

/*! @function
 @abstract  reset data to 0
 */
static void clear_counts(counts* p) {
    p->na = p->nc = p->ng = p->nn = p->nr = p->nt = p->nv = p->total = p->nx = 0;
    p->ref_b = 0;
    clear_str2int_hashmap(p->var);
    clear_str_set(p->umi);
}

/*! @function
 @abstract  reset data to 0
 */
static void clear_pcounts(pcounts* p) {
    clear_counts(p->mc);
    clear_counts(p->pc);
    p->pos = 0;
}

/*! @function
 @abstract convert variant hashmap to string for output "AC,AG,AT"
 */
static void get_var_string(str2intmap_t* vhash, char* out) {
    const char* reg;
    int z = 1;
    khint_t k;

    if (kh_size(*vhash) > 0) {
        for (k = kh_begin(*vhash); k < kh_end(*vhash); k++) {
            if (kh_exist(*vhash,k)) {
                reg = kh_key(*vhash,k);

                if (z == 1) {
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

/*! @function
 @abstract Add to count data
 @param  pd     processed pileup data structure
 @param  fi     file index
 @param  p      pointer to counts that are to be added to pd
 @param  ctig   contig
 @param  gpos   genomic base position
 @param  ref    genomic reference base
 @param  strand strand of base

 @return
 -1 on error, otherwise 0 if successful
 */
static int add_counts(PLP_DATA pd, int fi, counts* p, const char* ctig, int gpos, int ref, int strand) {

    int ret = 0;
    if (fi >= pd->nfiles || fi < 0) {
        REprintf("[raer internal] issue with file index %d", fi);
        return -1;
    }
    SEXP r = get_or_grow_PLP_DATA(pd, -1, PLP_DATA_LST), s;
    int li;
    int idx = pd->icnt;
    char* buf;
    char vout[12];
    get_var_string(&(p->var), vout);

    for (li = 0; li < LENGTH(r); ++li) {
        if (R_NilValue == (s = VECTOR_ELT(r, li)))
            continue;
        switch (li) {
        case SEQNAME_IDX:
            buf = R_Calloc(strlen(ctig) + 1, char);
            if (!buf) {
                REprintf("[raer internal ] add_plp_data: failed to allocate memory");
                return -1;
            }
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
            REprintf("[raer internal ] unhandled add_counts");
        ret = -1;
        break;
        }
    }
    return ret;
}

/*! @function
 @abstract Add to rowData stats
 @param  pd     processed pileup data structure
 @param  idx    vector index
 @param  rpbz   read position bias z-score
 @param  vdb    variant distance bias
 @param  sor    strand bias

 */
static void add_stats(PLP_DATA pd, int idx, double rpbz, double vdb, double sor) {
    // check if size is sufficient
    get_or_grow_PLP_DATA(pd, -1, SITE_DATA_LST);
    pd->sdat->rpbz[idx] = rpbz;
    pd->sdat->vdb[idx] = vdb;
    pd->sdat->sor[idx] = sor;
}


/*! @function
 @abstract Calculate site biases
 @param  pall   counts
 @param  res    vector containing statistics (length 6)

 @return
 -1 on error, 0 otherwise
 */
static int calc_biases(pall_counts* pall, double* res) {
    if (NBASE_POS != 100) return -1;
    if (pall->p_has_var) {
        res[0] = calc_mwu_biasZ(pall->p_ref_pos, pall->p_alt_pos, NBASE_POS, 0, 1);
        res[1] = calc_vdb(pall->p_alt_pos, NBASE_POS);
        res[2] = calc_sor(pall->s->ref_fwd, pall->s->ref_rev, pall->s->alt_fwd, pall->s->alt_rev);
    }
    if (pall->m_has_var) {
        res[3] = calc_mwu_biasZ(pall->m_ref_pos, pall->m_alt_pos, NBASE_POS, 0, 1);
        res[4] = calc_vdb(pall->m_alt_pos, NBASE_POS);
        res[5] = calc_sor(pall->s->ref_fwd, pall->s->ref_rev, pall->s->alt_fwd, pall->s->alt_rev);
    }
    return 0;
}

/*! @function
 @abstract Calculate site biases
 @param  pall   counts
 @param  res    vector containing statistics (length 6)

 @return
 -1 on error, 0 otherwise
 */
static int store_counts(PLP_DATA pd, pcounts* pc, const char* ctig,
                        const int gpos, int pref, int mref, double* stats,
                        mplp_conf_t* conf) {

    int i, nv, ret = 0;
    double vf;
    int p_has_v = 0;
    int m_has_v = 0;
    int p_is_ma = 0;
    int m_is_ma = 0;
    int write_p = 0;
    int write_m = 0;
    int write_only_v = 0;
    int vs = 0;
    khiter_t k;

    for (i = 0; i < conf->nbam; ++i) {
        if (conf->only_keep_variants[i]) {
            write_only_v = 1;
        }
    }

    // store counts if any samples are > min_depth and any samples have a variant
    // predicated on the true or false values in only_variants
    for (i = 0; i < conf->nbam; ++i) {
        // check depth
        if ((pc + i)->pc->total >= conf->min_depth && (pc + i)->pc->nv >= conf->ef.min_var_reads) {
            write_p = 1;
        }
        if ((pc + i)->mc->total >= conf->min_depth && (pc + i)->mc->nv >= conf->ef.min_var_reads) {
            write_m = 1;
        }

        // check allele freq, remove from hashmap if does not pass
        if (conf->min_af > 0) {
            if (kh_size((pc + i)->pc->var) > 0) {
                for (k = kh_begin((pc + i)->pc->var); k < kh_end((pc + i)->pc->var); k++) {
                    if (kh_exist((pc + i)->pc->var, k)) {
                        nv = kh_value((pc + i)->pc->var, k);
                        vf = (double) nv / (pc + i)->pc->total;
                        if (vf < conf->min_af) {
                            kh_del(str2intmap, (pc + i)->pc->var, k);
                        }
                    }
                }
            }

            if (kh_size((pc + i)->mc->var) > 0) {
                for (k = kh_begin((pc + i)->mc->var); k < kh_end((pc + i)->mc->var); k++) {
                    if (kh_exist((pc + i)->mc->var, k)) {
                        nv = kh_value((pc + i)->mc->var, k);
                        vf = (double) nv / (pc + i)->mc->total;
                        if (vf < conf->min_af) {
                            kh_del(str2intmap, (pc + i)->mc->var, k);
                        }
                    }
                }
            }
        }

        vs = kh_size((pc + i)->pc->var);
        if (vs > 0) {
            if (conf->only_keep_variants[i]) p_has_v = 1;
            if (vs > 1) p_is_ma = 1;
        }

        vs = kh_size((pc + i)->mc->var);
        if (vs > 0) {
            if (conf->only_keep_variants[i]) m_has_v = 1;
            if (vs > 1) m_is_ma = 1;
        }
    }

    if (write_only_v) {
        write_p = p_has_v && write_p;
        write_m = m_has_v && write_m;
    }

    if (!conf->report_multiallelics) {
        write_p = !p_is_ma && write_p;
        write_m = !m_is_ma && write_m;
    }

    if (write_p) {
        for (i = 0; i < conf->nbam; ++i) {
            ret = add_counts(pd, i, (pc + i)->pc, ctig, gpos, pref, '+');
        }
        add_stats(pd, pd->icnt, stats[0], stats[1], stats[2]);
        pd->icnt += 1;
    }

    if (write_m) {
        for (i = 0; i < conf->nbam; ++i) {
            ret = add_counts(pd, i, (pc + i)->mc, ctig, gpos, mref, '-');
        }
        add_stats(pd, pd->icnt, stats[3], stats[4], stats[5]);
        pd->icnt += 1;
    }
    return ret;
}

static int count_one_record(const bam_pileup1_t* p, pcounts* pc, mplp_conf_t* conf,
                            pall_counts* pall,  const char* ctig, int pos, int pref_b,
                            int mref_b, int read_base, int invert, int i) {
    int hret = 0, rpos = 0;
    char mvar[2];
    char pvar[2];
    mvar[1] = '\0';
    pvar[1] = '\0';
    khiter_t k;

    rpos = get_relative_position(p, NBASE_POS);
    if (rpos < 0) return -1;
    if (invert) {
        pc->mc->total += 1;

        // check if read base == ref base
        if (mref_b == read_base) {
            pc->mc->nr += 1;
            pall->m_ref_pos[rpos] += 1;
        } else {
            mvar[0] = read_base;

            // store variants as single bases in hash set
            char* var = strdup(mvar);

            k = kh_put(str2intmap, pc->mc->var, var, &hret);
            if (hret == 0) {
                free(var);
                kh_value(pc->mc->var, k) += 1;
            } else {
                kh_value(pc->mc->var, k) = 1;
            }

            pc->mc->nv += 1;
            pall->m_alt_pos[rpos] += 1;
            pall->m_has_var = 1;
        }
        switch (read_base) {
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

        if (pref_b == read_base) {
            pc->pc->nr += 1;
            pall->p_ref_pos[rpos] += 1;
        } else {

            pvar[0] = read_base;

            char* var = strdup(pvar);
            k = kh_put(str2intmap, pc->pc->var, var, &hret);
            if (hret == 0) {
                free(var);
                kh_value(pc->pc->var, k) += 1;
            } else {
                kh_value(pc->pc->var, k) = 1;
            }

            pc->pc->nv += 1;
            pall->p_alt_pos[rpos] += 1;
            pall->p_has_var = 1;
        }
        switch (read_base) {
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
static int check_read_filters(const bam_pileup1_t* p, mplp_conf_t* conf, int baq, int maq) {
    int res = 0;
    // skip indel and ref skip ;
    if (p->is_del || p->is_refskip) return (2) ;

    // filter based on mapq as not able to filter per file in pileup
    if (p->b->core.qual < maq) return (1) ;

    // check base quality
    int bq = p->qpos < p->b->core.l_qseq
    ? bam_get_qual(p->b)[p->qpos]
    : 0;
    if (bq < baq) {
        if (bq == 0) { // overlapping read pairs get base qualities set to 0 by mplp
            return (2);
        } else {
            return (1);
        }
    }

    // check if pos is within x dist from 5' end of read, qpos is 0-based
    if (conf->ef.trim_f5p > 0 || conf->ef.trim_f3p > 0) {
        res = check_variant_fpos(p->b, p->qpos, conf->ef.trim_f5p, conf->ef.trim_f3p);
        if (res < 0) return -1; // error
        if (res > 0) return 1;
    }

    if (conf->ef.trim_i5p > 0 || conf->ef.trim_i3p > 0) {
        res = check_variant_pos(p->b, p->qpos, conf->ef.trim_i5p, conf->ef.trim_i3p);
        if (res < 0) return -1; // error
        if (res > 0) return 1;
    }

    // check for splice in alignment nearby
    if (conf->ef.splice_dist && dist_to_splice(p->b, p->qpos, conf->ef.splice_dist) >= 0) return (1);

    // check if site in splice overhang and > min_overhang
    if (conf->ef.min_overhang) {
        res = check_splice_overhang(p->b, p->qpos, conf->ef.min_overhang);
        if (res == -2) return -1; // error
        if (res > 0) return (1);
    }

    // check if indel event nearby
    if (conf->ef.indel_dist && dist_to_indel(p->b, p->qpos, conf->ef.indel_dist) >= 0) return (1);

    return 0;
}

/* return codes,
 * -1 = read has no UMI,
 * 0 = read has UMI already seen
 * 1 = read has new UMI
 */
static int check_umi(const bam_pileup1_t* p, mplp_conf_t* conf,
                     pcounts* plpc, int invert) {
    char* umi;
    int ret;
    if (invert) {
        umi = get_aux_ztag(p->b, conf->umi_tag);
        if (umi == NULL) return -1;
        char* umi_val = strdup(umi);
        kh_put(strset, plpc->mc->umi, umi_val, &ret);
        if (ret == 0) {
            free(umi_val);
            return (0);
        }
    } else {
        umi = get_aux_ztag(p->b, conf->umi_tag);
        if (umi == NULL) return -1;
        char* umi_val = strdup(umi);
        kh_put(strset, plpc->pc->umi, umi_val, &ret);
        if (ret == 0) {
            free(umi_val);
            return (0);
        }
    }
    return (1);
}

static int run_pileup(char** cbampaths, char** cindexes,
                      PLP_DATA pd, mplp_conf_t* conf) {

    hts_set_log_level(HTS_LOG_ERROR);

    mplp_aux_t** data;
    int i, tid, *n_plp, ret = 0;
    hts_pos_t pos, beg0 = 0, end0 = INT32_MAX, ref_len;

    const bam_pileup1_t** plp;
    mplp_ref_t mp_ref = MPLP_REF_INIT;
    bam_mplp_t iter;
    bam_hdr_t* h = NULL; /* header of first file in input list */
    char* ref;

    data = calloc(conf->nbam, sizeof(mplp_aux_t*));
    plp = calloc(conf->nbam, sizeof(bam_pileup1_t*));
    n_plp = calloc(conf->nbam, sizeof(int));

    pcounts* plpc;
    plpc = R_Calloc(conf->nbam, pcounts);
    for (i = 0; i < conf->nbam; ++i) {
        plpc[i].pc = R_Calloc(1, counts);
        plpc[i].mc = R_Calloc(1, counts);

        //initialize variant string set (AG, AT, TC) etc.
        if (plpc[i].pc->var == NULL) plpc[i].pc->var = kh_init(str2intmap);
        if (plpc[i].mc->var == NULL) plpc[i].mc->var = kh_init(str2intmap);

        //initialize umi string set
        if (plpc[i].pc->umi == NULL) plpc[i].pc->umi = kh_init(strset);
        if (plpc[i].mc->umi == NULL) plpc[i].mc->umi = kh_init(strset);
    }

    pall_counts* pall = init_pall_counts();

    double* site_stats = R_Calloc(6, double);

    // read the header of each file in the list and initialize data
    for (i = 0; i < conf->nbam; ++i) {
        bam_hdr_t* h_tmp;
        data[i] = calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open(cbampaths[i], "rb");
        if (!data[i]->fp) {
            Rf_error("[raer internal] failed to open %s: %s\n",
                     cbampaths[i], strerror(errno));
        }
        if (conf->fai_fname) {
            if (hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
                Rf_error("[raer internal] failed to process %s: %s\n",
                         conf->fai_fname, strerror(errno));
            }
        }
        data[i]->conf = conf;
        data[i]->ref = &mp_ref;

        h_tmp = sam_hdr_read(data[i]->fp);
        if (!h_tmp) {
            Rf_error("[raer internal] fail to read the header of %s\n", cbampaths[i]);
        }

        if (conf->reg) {
            hts_idx_t* idx = NULL;
            idx = sam_index_load2(data[i]->fp, cbampaths[i], cindexes[i]) ;

            if (idx == NULL) {
                Rf_error("[raer internal] fail to load bamfile index for %s\n", cbampaths[i]);
            }

            if ((data[i]->iter=sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
                Rf_error("[raer internal] fail to parse region '%s' with %s\n", conf->reg, cbampaths[i]);
            }

            if (i == 0) {
                beg0 = data[i]->iter->beg;
                end0 = data[i]->iter->end;
            }
            hts_idx_destroy(idx);
        } else {
            data[i]->iter = NULL;
        }

        if (i == 0) {
            h = data[i]->h = h_tmp;
        } else {
            bam_hdr_destroy(h_tmp);
            data[i]->h = h;
        }
    }

    iter = bam_mplp_init(conf->nbam, readaln, (void**)data);

    ret = bam_mplp_init_overlaps(iter) ;
    if (ret < 0) {
        REprintf("[raer internal] issue initializing iterator");
        ret = -1;
        goto fail;
    }

    // set max depth
    bam_mplp_set_maxcnt(iter, conf->max_depth);

    if (!iter) {
        REprintf("[raer internal] issue setting max depth on iterator");
        ret = -1;
        goto fail;
    }


    int n_iter = 0;
    while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
        ++n_iter;

        if (conf->reg && (pos < beg0 || pos >= end0)) continue;

        mplp_get_ref(data[0], tid, &ref, &ref_len);

        if (tid < 0) break;

        // check user interrupt, using a 2^k value is faster
        if (n_iter % 262144 == 0) {
            if(checkInterrupt()){
                REprintf("[raer internal] user interrupt detected, exiting\n");
                ret = -1;
                goto fail;
            }
        }

        // ensure position is in requested intervals
        if (conf->reg_idx && tid >= 0) {
            int ol = regidx_overlap(conf->reg_idx,
                                    sam_hdr_tid2name(h, tid),
                                    pos,
                                    pos,
                                    conf->reg_itr);
            if (!ol) continue;
        }

        // get reference base on +/- strand
        int pref_b, mref_b;
        pref_b = (ref && pos < ref_len)? ref[pos] : 'N' ;
        pref_b = toupper(pref_b);

        if (pref_b == 'N') continue;

        mref_b = comp_base[(unsigned char) pref_b];

        // reset count structure
        clear_pall_counts(pall);
        memset(site_stats, 0, sizeof(double) * 6);
        for (i = 0; i < conf->nbam; ++i) {
            clear_pcounts(&plpc[i]);
        }

        // check if site is in a homopolymer
        if (conf->ef.nmer > 0 && check_simple_repeat(&ref, &ref_len, pos, conf->ef.nmer)) continue;

        // check if read count less than min_depth
        int pass_reads = 0;
        for (i = 0; i < conf->nbam; ++i) {
            if (n_plp[i] >= conf->min_depth) {
                pass_reads = 1;
            }
        }
        if (!pass_reads) continue;

        // iterate through bam files
        for (i = 0; i < conf->nbam; ++i) {
            int j;
            // iterate through reads that overlap position
            for (j = 0; j < n_plp[i]; ++j) {

                const bam_pileup1_t* p = plp[i] + j;

                // get read base
                int c = p->qpos < p->b->core.l_qseq
                ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
                : 'N';
                if (c == 'N') continue;

                int invert = invert_read_orientation(p->b, conf->libtype[i]);
                if (invert < 0) {
                    REprintf("[raer internal] invert read orientation failure %i\n", invert);
                    ret = -1;
                    goto fail;
                }

                // remove bad reads
                int rret = check_read_filters(p, conf, conf->min_bq, conf->min_mqs[i]);

                // only keep first read with a UMI tag per position
                if (conf->umi) {
                    int uret = check_umi(p, conf, &plpc[i], invert);
                    if (uret != 1) continue;
                }

                if (rret > 0) {
                    if (rret == 1) {
                        if (invert) {
                            plpc[i].mc->nx += 1;
                        } else {
                            plpc[i].pc->nx += 1;
                        }
                    }
                    continue;
                }
                int ci = c;
                if (invert) ci = (char)comp_base[(unsigned char)c];

                // check read for >= mismatch different types and at least n_mm mismatches
                if (conf->ef.n_mm_type > 0 || conf->ef.n_mm > 0) {
                    if ((invert && mref_b != ci) || pref_b != ci) {
                        int m = parse_mismatches(p->b, conf->ef.n_mm_type, conf->ef.n_mm);
                        if (m == -1) {
                            ret = -1;
                            goto fail;
                        } else if (m > 0) {
                            continue;
                        }
                    }
                }


                if (pref_b == c) {
                    if (bam_is_rev(p->b)) {
                        pall->s->ref_rev += 1;
                    } else {
                        pall->s->ref_fwd += 1;
                    }
                } else {
                    if (bam_is_rev(p->b)) {
                        pall->s->alt_rev += 1;
                    } else {
                        pall->s->alt_fwd += 1;
                    }
                }

                // increment counts per sample
                int cret = count_one_record(p, &plpc[i], conf, pall,
                                            sam_hdr_tid2name(h, tid), pos,
                                            pref_b, mref_b,
                                            ci, invert, i);
                if (cret < 0) {
                    ret = -1;
                    goto fail;
                }

            }
        }

        int sres = calc_biases(pall, site_stats);
        if (sres < 0) {
            ret = -1;
            goto fail;
        }

        // write or store records if pass depth criteria
        sres = store_counts(pd, plpc, sam_hdr_tid2name(h, tid),
                            pos, pref_b, mref_b, site_stats, conf);
        if (sres < 0) {
            REprintf("[raer internal] failed storing counts, %s %d\n",
                     sam_hdr_tid2name(h, tid),
                     pos);
            ret = -1;
            goto fail;
        }
    }

    fail:
        if (ret >= 0) finish_PLP_DATA(pd);

        bam_mplp_destroy(iter);
        bam_hdr_destroy(h);

        for (i = 0; i < conf->nbam; ++i) {
            sam_close(data[i]->fp);
            if (data[i]->iter) hts_itr_destroy(data[i]->iter);
            free(data[i]);
            clear_pcounts(&plpc[i]);
            kh_destroy(str2intmap, plpc[i].mc->var);
            kh_destroy(str2intmap, plpc[i].pc->var);
            kh_destroy(strset, plpc[i].mc->umi);
            kh_destroy(strset, plpc[i].pc->umi);
            R_Free(plpc[i].pc);
            R_Free(plpc[i].mc);

        }

        free(data);
        free(plp);
        free(n_plp);
        free(mp_ref.ref[0]);
        free(mp_ref.ref[1]);

        if (pall) {
            R_Free(pall->p_ref_pos);
            R_Free(pall->p_alt_pos);
            R_Free(pall->m_alt_pos);
            R_Free(pall->m_alt_pos);
            R_Free(pall);
        }

        R_Free(plpc);
        if (pd->fps) R_Free(pd->fps);
        R_Free(pd->pdat);

        if (site_stats) R_Free(site_stats);

    return ret;
}


static void check_plp_args(SEXP bampaths, SEXP indexes, SEXP n, SEXP fapath, SEXP region,
                           SEXP int_args, SEXP dbl_args, SEXP lgl_args,
                           SEXP libtype, SEXP only_keep_variants, SEXP min_mapQ,
                           SEXP umi) {
    if (!IS_INTEGER(n) || (LENGTH(n) != 1)) {
        Rf_error("'n' must be integer(1)");
    }
    int n_files = INTEGER(n)[0];

    if (!IS_CHARACTER(bampaths) || (LENGTH(bampaths) != n_files)) {
        Rf_error("'bampaths' must be character vector equal in length to number of bam files");
    }

    if (!IS_CHARACTER(indexes) || (LENGTH(indexes) != n_files)) {
        Rf_error("'indexes' must be character vector equal in length to number of bam files");
    }

    if (!IS_CHARACTER(fapath) || (LENGTH(fapath) != 1)) {
        Rf_error("'fapath' must be character(1)");
    }

    if (!IS_CHARACTER(region) || (LENGTH(region) > 1)) {
        Rf_error("'region' must be character of length 0 or 1");
    }

    // vectors populated with parameters of fixed sizes
    if (!IS_INTEGER(int_args) || (LENGTH(int_args) != 14)) {
        Rf_error("'int_args' must be integer of length 14");
    }

    if (!IS_NUMERIC(dbl_args) || (LENGTH(dbl_args) != 5)) {
        Rf_error("'dbl_args' must be numeric of length 5");
    }

    if (!IS_LOGICAL(lgl_args) || (LENGTH(lgl_args) != 1)) {
        Rf_error("'lgl_args' must be logical of length 1");
    }

    // args depending on n_files
    if (!IS_INTEGER(libtype) || (LENGTH(libtype) != n_files)) {
        Rf_error("'lib_type' must be integer of same length as bamfiles");
    }

    if (!IS_LOGICAL(only_keep_variants) || (LENGTH(only_keep_variants) != n_files)) {
        Rf_error("'only_keep_variants' must be logical of same length as bamfiles");
    }

    if (!IS_INTEGER(min_mapQ) || (LENGTH(min_mapQ) != n_files)) {
        Rf_error("'min_mapQ' must be integer of same length as bamfiles");
    }

    // other args

    if (!IS_CHARACTER(umi) || (LENGTH(umi) > 1)) {
        Rf_error("'umi' must be character of length 0 or 1");
    }

}

static int set_mplp_conf(mplp_conf_t* conf, int n_bams,
                         char* fafn, char* qregion, regidx_t* idx,
                         int* i_args, double* d_args, int* b_args,
                         int* libtypes, int* keep_variants, int* min_mapqs,
                         char* umi) {
    int ret = 0;
    if (n_bams <= 0) {
        REprintf("[raer internal] invalid bam input");
        return -1;
    }
    conf->nbam = n_bams;
    conf->nfps = n_bams + 1; // 1 extra file stores the rowdata information

    if (fafn) {
        conf->fai_fname = fafn;
        conf->fai = fai_load(conf->fai_fname);
        if (conf->fai == NULL) {
            REprintf("[raer internal] unable to load fasta index");
            return -1;
        }
    }

    // single region for pileup
    if (qregion) conf->reg = qregion;

    // store interval index
    if (idx) {
        conf->reg_idx = (regidx_t*)idx;
        conf->reg_itr = regitr_init(conf->reg_idx);
        if (!conf->reg_itr) {
            REprintf("[raer internal] unable to load region index iterator");
            return -1;
        }
    }

    memset(&conf->ef, 0, sizeof(efilter_t));

    conf->max_depth        = i_args[0];
    conf->min_depth        = i_args[1];
    conf->min_bq           = i_args[2];
    conf->ef.trim_i5p      = i_args[3];
    conf->ef.trim_i3p      = i_args[4];
    conf->ef.indel_dist    = i_args[5];
    conf->ef.splice_dist   = i_args[6];
    conf->ef.min_overhang  = i_args[7];
    conf->ef.nmer          = i_args[8];
    conf->ef.min_var_reads = i_args[9];
    conf->ef.n_mm_type     = i_args[10];
    conf->ef.n_mm          = i_args[11];
    conf->keep_flag[0]     = i_args[12];
    conf->keep_flag[1]     = i_args[13];

    conf->ef.trim_f5p      = d_args[0];
    conf->ef.trim_f3p      = d_args[1];
    conf->min_af           = d_args[2];
    conf->read_qual.pct    = d_args[3];
    conf->read_qual.minq   = d_args[4];

    conf->report_multiallelics = b_args[0];

    conf->libtype            = libtypes;
    conf->only_keep_variants = keep_variants;

    // if multiple bam files, use minimum mapQ for initial filtering,
    // then later filter in pileup loop
    conf->min_mqs = min_mapqs;
    if (min_mapqs[0] > 0) {
        int mq;
        mq = min_mapqs[0];
        for (int i = 0; i < n_bams; ++i) {
            if (min_mapqs[i] < mq) {
                mq = min_mapqs[i];
            }
        }
        conf->min_global_mq = mq;
    }

    if (umi) {
        conf->umi = 1;
        conf->umi_tag = umi;
    }

    return ret;
}

SEXP pileup(SEXP bampaths, SEXP indexes, SEXP n, SEXP fapath, SEXP region, SEXP lst,
            SEXP int_args, SEXP dbl_args, SEXP lgl_args,
            SEXP libtype, SEXP only_keep_variants, SEXP min_mapQ, SEXP umi) {

    check_plp_args(bampaths, indexes, n, fapath, region,
                   int_args, dbl_args, lgl_args,
                   libtype, only_keep_variants, min_mapQ, umi);

    regidx_t* idx = NULL;
    if (!Rf_isNull(lst) && Rf_length(VECTOR_ELT(lst, 0)) > 0) {
        idx = regidx_build(lst, 2);
        if (!idx) Rf_error("Failed to build region index");
    }

    int i, ret = 0;;
    char** cbampaths = (char**) R_alloc(sizeof(const char*), Rf_length(bampaths));
    char** cindexes = (char**) R_alloc(sizeof(const char*), Rf_length(indexes));
    for (i = 0; i < LENGTH(bampaths); ++i) {
        cbampaths[i] = (char*) translateChar(STRING_ELT(bampaths, i));
        cindexes[i] = (char*) translateChar(STRING_ELT(indexes, i));
    }

    char* cfafn = (char*) translateChar(STRING_ELT(fapath, 0));

    char* cregion = LENGTH(region) == 0 ?
    NULL : (char*) translateChar(STRING_ELT(region, 0));

    char* umi_tag = LENGTH(umi) == 0 ?
    NULL : (char*) translateChar(STRING_ELT(umi, 0));

    mplp_conf_t ga;
    memset(&ga, 0, sizeof(mplp_conf_t));

    int nbam = INTEGER(n)[0];
    ret = set_mplp_conf(&ga, nbam, cfafn, cregion, idx,
                        INTEGER(int_args), REAL(dbl_args), LOGICAL(lgl_args),
                        INTEGER(libtype), LOGICAL(only_keep_variants), INTEGER(min_mapQ),
                        umi_tag);
    SEXP result;
    result = PROTECT(pileup_result_init(nbam));
    PLP_DATA pd;
    if (ret >= 0) {
        pd = init_PLP_DATA(result,  nbam);
        ret = run_pileup(cbampaths, cindexes, pd, &ga);
    }

    // clean up
    if (ga.fai) fai_destroy(ga.fai);
    if (ga.reg_itr) regitr_destroy(ga.reg_itr);
    if (ga.reg_idx) regidx_destroy(ga.reg_idx);

    UNPROTECT(1);
    if (ret < 0) Rf_error("[raer internal] error detected during pileup");

    return result ;
}
