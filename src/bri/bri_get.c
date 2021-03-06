//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// bri - simple utility to provide random access to
//       bam records by read name
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "bri_index.h"


#include <Rinternals.h>
//
// Getopt
//
// enum {
//     OPT_HELP = 1,
// };
//
// static const char* shortopts = ":i:"; // placeholder
// static const struct option longopts[] = {
//     { "help",                      no_argument,       NULL, OPT_HELP },
//     { "index",               required_argument,       NULL,      'i' },
//     { NULL, 0, NULL, 0 }
// };

void print_usage_get()
{
    Rf_error("usage: bri get [-i <index_filename.bri>] <input.bam> <output.bam> <readname>\n");
}

// comparator used by bsearch, direct strcmp through the name pointer
int compare_records_by_readname_ptr(const void* r1, const void* r2)
{
    const char* n1 = ((bam_read_idx_record*)r1)->read_name.ptr;
    const char* n2 = ((bam_read_idx_record*)r2)->read_name.ptr;
    return strcmp(n1, n2);
}

//
void bam_read_idx_get_range(const bam_read_idx* bri, const char* readname, bam_read_idx_record** start, bam_read_idx_record** end)
{
    // construct a query record to pass to bsearch
    bam_read_idx_record query;
    query.read_name.ptr = readname;
    query.file_offset = 0;

    // if rec is NULL then readname does not appear in index
    bam_read_idx_record* rec =
        bsearch(&query, bri->records, bri->record_count, sizeof(bam_read_idx_record), compare_records_by_readname_ptr);
    if(rec == NULL) {
        *start = NULL;
        *end = NULL;
        return;
    }

    // rec points to a valid record, but it can be an arbitrary record in the range
    // move start to the first record in the range, and end to be one past the end
    size_t sri, eri;
    sri = eri = rec - bri->records;
    assert(bri->records[sri].file_offset == rec->file_offset);

    while(sri > 0 && bri->records[sri].read_name.ptr == bri->records[sri - 1].read_name.ptr) {
        sri -= 1;
    }
    assert(strcmp(bri->records[sri].read_name.ptr, readname) == 0);

    do {
        eri += 1;
    } while(eri < bri->record_count && bri->records[eri].read_name.ptr == bri->records[sri].read_name.ptr);
    assert(eri == bri->record_count || strcmp(bri->records[eri].read_name.ptr, readname) != 0);
    //Rf_error("r: %zu sri: %zu eri: %zu\n", rec - bri->records, sri, eri);
    *start = &bri->records[sri];
    *end = &bri->records[eri];
}

//
void bam_read_idx_get_by_record(htsFile* fp, bam_hdr_t* hdr, bam1_t* b, bam_read_idx_record* bri_record)
{
    int ret = bgzf_seek(fp->fp.bgzf, bri_record->file_offset, SEEK_SET);
    if(ret != 0) {
        Rf_error("[bri] bgzf_seek failed\n");
    }

    ret = sam_read1(fp, hdr, b);
    if(ret < 0) {
        Rf_error("[bri] sam_read1 failed\n");
    }
}


int bam_read_idx_get(const char* input_bam,
                     const char* output_bam,
                     const char** bcs,
                     int nbcs) {
  char* input_bri = NULL;
  bam_read_idx* bri = bam_read_idx_load(input_bam, input_bri);

  htsFile *bam_fp = sam_open(input_bam, "r");
  bam_hdr_t *h = sam_hdr_read(bam_fp);
  htsFile *out_fp ;
  out_fp = hts_open(output_bam, "wb");

  int ret = 0;
  ret = sam_hdr_write(out_fp, h);
  if(ret != 0){
    Rf_error("[bri] sam_hdr_write failed\n");
  }

  bam_read_idx_record* start;
  bam_read_idx_record* end;

  for(int i = 0; i < nbcs; i++) {
    const char* bc = bcs[i];
    bam_read_idx_get_range(bri, bc, &start, &end);

    bam1_t *b = bam_init1();
    int n_rec = 0;
    while(start != end) {
      int ret = bgzf_seek(bam_fp->fp.bgzf , start->file_offset, SEEK_SET);
      if(ret != 0) {
        Rf_error("[bri] bgzf_seek failed\n");
      }

      while(n_rec < start->n_aln){
        ret = sam_read1(bam_fp, h, b);
        if(ret < 0) {
          Rf_error("[bri] sam_read1 failed\n");
        }
        int ret = sam_write1(out_fp, h, b);
        if(ret < 0) {
          Rf_error("[bri] sam_write1 failed\n");
        }
        n_rec += 1;

      }
      start++;
    }
    bam_destroy1(b);
  }

  hts_close(out_fp);
  bam_hdr_destroy(h);
  hts_close(bam_fp);
  bam_read_idx_destroy(bri);
  bri = NULL;

  return 0;


}

//
// int bam_read_idx_get_main(int argc, char** argv)
// {
//     char* input_bri = NULL;
//
//     int die = 0;
//     for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
//         switch (c) {
//             case OPT_HELP:
//                 print_usage_get();
//                 exit(EXIT_SUCCESS);
//             case 'i':
//                 input_bri = optarg;
//         }
//     }
//
//     if (argc - optind < 2) {
//         Rf_error("bri get: not enough arguments\n");
//         die = 1;
//     }
//
//     if(die) {
//         print_usage_get();
//     }
//
//     char* input_bam = argv[optind++];
//     char* output_bam = argv[optind++];
//     bam_read_idx* bri = bam_read_idx_load(input_bam, input_bri);
//
//     htsFile *bam_fp = sam_open(input_bam, "r");
//     bam_hdr_t *h = sam_hdr_read(bam_fp);
//     htsFile *out_fp ;
//     out_fp = hts_open(output_bam, "wb");
//
//     int ret = 0;
//     ret = sam_hdr_write(out_fp, h);
//     if(ret != 0){
//       Rf_error("[bri] sam_hdr_write failed\n");
//     }
//
//     bam_read_idx_record* start;
//     bam_read_idx_record* end;
//
//     for(int i = optind; i < argc; i++) {
//         char* readname = argv[i];
//         bam_read_idx_get_range(bri, readname, &start, &end);
//
//         bam1_t *b = bam_init1();
//         int n_rec = 0;
//         while(start != end) {
//             int ret = bgzf_seek(bam_fp->fp.bgzf , start->file_offset, SEEK_SET);
//             if(ret != 0) {
//                 Rf_error("[bri] bgzf_seek failed\n");
//             }
//
//             while(n_rec < start->n_aln){
//                 ret = sam_read1(bam_fp, h, b);
//                 if(ret < 0) {
//                     Rf_error("[bri] sam_read1 failed\n");
//
//                 }
//                 int ret = sam_write1(out_fp, h, b);
//                 if(ret < 0) {
//                     Rf_error("[bri] sam_write1 failed\n");
//                 }
//                 n_rec += 1;
//
//             }
//             start++;
//         }
//         bam_destroy1(b);
//     }
//
//     hts_close(out_fp);
//     bam_hdr_destroy(h);
//     hts_close(bam_fp);
//     bam_read_idx_destroy(bri);
//     bri = NULL;
//
//     return 0;
// }
