//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// bri - simple utility to provide random access to
//       bam records by read name
//

// avoid warnings in qsort_r
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "bri_index.h"
#include "sort_r.h"

#include <Rinternals.h>

//#define BRI_INDEX_DEBUG 1
char verbose = 0;
// make the index filename based on the name of input_bam
// caller must free the returned pointer
char* generate_index_filename(const char* input_bam, const char* input_bri)
{
    char* out_fn;

    if(input_bri != NULL) {
        out_fn = malloc(strlen(input_bri) + 1);
        if(out_fn == NULL) {
            Rf_error("error making index: %s", input_bri);
        }
        strcpy(out_fn, input_bri);
    } else {
        out_fn = malloc(strlen(input_bam) + 5);
        if(out_fn == NULL) {
          Rf_error("error making index: %s", input_bri);
        }
        strcpy(out_fn, input_bam);
        strcat(out_fn, ".bri");
    }

    return out_fn;
}

//
bam_read_idx* bam_read_idx_init()
{
    bam_read_idx* bri = (bam_read_idx*)malloc(sizeof(bam_read_idx));

    bri->name_capacity_bytes = 0;
    bri->name_count_bytes = 0;
    bri->readnames = NULL;

    bri->record_capacity = 0;
    bri->record_count = 0;
    bri->records = NULL;

    return bri;
}

//
void bam_read_idx_destroy(bam_read_idx* bri)
{
#ifdef BRI_INDEX_DEBUG
    REprintf("[bri-destroy] %zu name bytes %zu records\n", bri->name_count_bytes, bri->record_count);
#endif
    free(bri->readnames);
    bri->readnames = NULL;

    free(bri->records);
    bri->records = NULL;

    free(bri);
}

// comparison function used by quicksort, names pointers to a block of C strings
// and allows us to indirectly sorted records by their offset.
int compare_records_by_readname_offset(const void* r1, const void* r2, void* names)
{
    const char* cnames = (const char*)names;
    const char* n1 = cnames + ((bam_read_idx_record*)r1)->read_name.offset;
    const char* n2 = cnames + ((bam_read_idx_record*)r2)->read_name.offset;
    return strcmp(n1, n2);
}

//
void bam_read_idx_save(bam_read_idx* bri, const char* filename)
{
    FILE* fp = fopen(filename, "wb");
    if(!fp) {
      Rf_error("error writing index: %s", filename);
    }

    // Sort records by readname
    sort_r(bri->records, bri->record_count, sizeof(bam_read_idx_record), compare_records_by_readname_offset, bri->readnames);

    // write header, containing file version, the size (in bytes) of the read names
    // and the number of records. The readnames size is a placeholder and will be
    // corrected later.

    // version
    size_t FILE_VERSION = 1;
    fwrite(&FILE_VERSION, sizeof(FILE_VERSION), 1, fp);

    // readname length
    size_t readname_bytes = 0;
    fwrite(&readname_bytes, sizeof(readname_bytes), 1, fp);

    // num records
    fwrite(&bri->record_count, sizeof(bri->record_count), 1, fp);

    // Pass 1: count up the number of non-redundant read names, write them to disk
    // Also store the position in the file where the read name for each record was written
    size_t* disk_offsets_by_record = malloc(bri->record_count * sizeof(size_t));
    const char* rn = bri->readnames; // for convenience

    for(size_t i = 0; i < bri->record_count; ++i) {

        int redundant = i > 0 &&
            strcmp(rn + bri->records[i].read_name.offset, rn + bri->records[i - 1].read_name.offset) == 0;

        if(!redundant) {
            disk_offsets_by_record[i] = readname_bytes; // current position in file
            size_t len = strlen(rn + bri->records[i].read_name.offset) + 1;
            fwrite(rn + bri->records[i].read_name.offset, len, 1, fp);
            readname_bytes += len;
        } else {
            disk_offsets_by_record[i] = disk_offsets_by_record[i - 1];
        }

#ifdef BRI_INDEX_DEBUG
        REprintf("record %zu name: %s redundant: %d do: %zu offset: %zu\n",
            i, bri->readnames + bri->records[i].read_name.offset, redundant, disk_offsets_by_record[i], bri->records[i].file_offset);
#endif
    }

    // Pass 2: write the records, getting the read name offset from the disk offset (rather than
    // the memory offset stored)
    for(size_t i = 0; i < bri->record_count; ++i) {
        bam_read_idx_record brir = bri->records[i];
        brir.read_name.offset = disk_offsets_by_record[i];
        fwrite(&brir, sizeof(brir), 1, fp);
#ifdef BRI_INDEX_DEBUG
        REprintf("[bri-save] record %zu %s name offset: %zu file offset: %zu\n",
            i, bri->readnames + bri->records[i].read_name.offset, disk_offsets_by_record[i], bri->records[i].file_offset);
#endif
    }

    // finish by writing the actual size of the read name segment
    fseek(fp, sizeof(FILE_VERSION), SEEK_SET);
    fwrite(&readname_bytes, sizeof(readname_bytes), 1, fp);

    free(disk_offsets_by_record);
    fclose(fp);
}

// add a record to the index, growing the dynamic arrays as necessary
void bam_read_idx_add(bam_read_idx* bri, const char* readname, size_t offset)
{
    //
    // add readname to collection if not seen, otherwise keep original offset
    //

    if(bri->record_count > 0){
        const char* orn = bri->readnames + bri->records[bri->record_count - 1].read_name.offset ;
        if(strcmp(orn, readname) == 0) {
            bri->records[bri->record_count - 1].n_aln += 1;
          //  bri->records[bri->record_count - 1].file_offset = offset;
           // printf("%s %s\n", bri->readnames + bri->records[bri->record_count - 1].read_name.offset, readname);
            return;
        }
    }

    size_t len = strlen(readname) + 1;
    if(bri->name_capacity_bytes <= bri->name_count_bytes + len) {

        // if already allocated double size, if initialization start with 1Mb
        bri->name_capacity_bytes = bri->name_capacity_bytes > 0 ? 2 * bri->name_capacity_bytes : 1024*1024;
#ifdef BRI_INDEX_DEBUG
        REprintf("[bri] allocating %zu bytes for names\n", bri->name_capacity_bytes);
#endif

        bri->readnames = realloc(bri->readnames, bri->name_capacity_bytes);
        if(bri->readnames == NULL) {
            Rf_error("[bri] malloc failed\n");
        }
    }

    // in principle this incoming name can be so larger than doubling the size
    // doesn't allow it to fit. This really shouldn't happen so we'll just exit here
    if(bri->name_capacity_bytes <= bri->name_count_bytes + len) {
        Rf_error("[bri] incoming name with length %zu is too large (%zu %zu)\n", len, bri->name_count_bytes + len, bri->name_capacity_bytes);
    }

    // copy name
    size_t name_offset = bri->name_count_bytes;
    strncpy(bri->readnames + bri->name_count_bytes, readname, len);
    bri->name_count_bytes += len;
    assert(bri->readnames[bri->name_count_bytes - 1] == '\0');

    //
    // add record
    //
    if(bri->record_count == bri->record_capacity) {
        bri->record_capacity = bri->record_capacity > 0 ? 2 * bri->record_capacity : 1024;
        bri->records = realloc(bri->records, sizeof(bam_read_idx_record) * bri->record_capacity);
        if(bri->records == NULL) {
            Rf_error("[bri] malloc failed\n");
        }
    }
    bri->records[bri->record_count].n_aln = 1;
    bri->records[bri->record_count].read_name.offset = name_offset;
    bri->records[bri->record_count].file_offset = offset;
    bri->record_count += 1;

}

//
void bam_read_idx_build(const char* filename,
                        const char* output_bri,
                        const char* tag)
{
    htsFile *fp = hts_open(filename, "r");
    if(fp == NULL) {
        Rf_error("[bri] could not open %s\n", filename);
    }

    bam_read_idx* bri = bam_read_idx_init();

    bam1_t* b = bam_init1();
    bam_hdr_t *h = sam_hdr_read(fp);
    int ret = 0;
    int niter = 0;
    size_t file_offset = bgzf_tell(fp->fp.bgzf);
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        char* readname ;
        uint8_t* aux_info = bam_aux_get(b, tag) ;
        if (aux_info) {
            readname = bam_aux2Z(aux_info) ;
        } else {
            continue;
        }

        bam_read_idx_add(bri, readname, file_offset);

        bam_read_idx_record brir = bri->records[bri->record_count - 1];
        if(verbose && (niter == 1 || niter % 100000 == 0)) {
            REprintf("[bri-build] record %zu [%zu %zu] %s\n",
                bri->record_count,
                brir.read_name.offset,
                brir.file_offset,
                bri->readnames + brir.read_name.offset
            );
        }
        niter += 1;
        // update offset for next record
        file_offset = bgzf_tell(fp->fp.bgzf);
    }

    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);

    // save to disk and cleanup
    if(verbose) {
        REprintf("[bri-build] writing to disk...\n");
    }

    char* out_fn = generate_index_filename(filename, output_bri);
    bam_read_idx_save(bri, out_fn);

    if(verbose) {
        REprintf("[bri-build] wrote index for %zu records.\n", bri->record_count);
    }

    free(out_fn);
    bam_read_idx_destroy(bri);
}

void print_error_and_exit(const char* msg)
{
    Rf_error("[bri] %s\n", msg);
}

//
bam_read_idx* bam_read_idx_load(const char* input_bam, const char* input_bri)
{
    char* index_fn = generate_index_filename(input_bam, input_bri);
    FILE* fp = fopen(index_fn, "rb");
    if(fp == NULL) {
        Rf_error("[bri] index file not found for %s\n", input_bam);
    }

    bam_read_idx* bri = bam_read_idx_init();
    size_t file_version;
    // currently ignored
    int bytes_read = fread(&file_version, sizeof(file_version), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // read size of readames segment
    bytes_read = fread(&bri->name_count_bytes, sizeof(bri->name_count_bytes), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }
    bri->name_capacity_bytes = bri->name_count_bytes;

    // read number of records on disk
    bytes_read = fread(&bri->record_count, sizeof(bri->record_count), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    bri->record_capacity = bri->record_count;

    // allocate filenames
    bri->readnames = malloc(bri->name_capacity_bytes);
    if(bri->readnames == NULL) {
        Rf_error("[bri] failed to allocate %zu bytes for read names\n", bri->name_capacity_bytes);
    }

    // allocate records
    bri->records = malloc(bri->record_capacity * sizeof(bam_read_idx_record));

    // read the names
    bytes_read = fread(bri->readnames, bri->name_count_bytes, 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // read the records
    bytes_read = fread(bri->records, bri->record_count, sizeof(bam_read_idx_record), fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // convert read name offsets to direct pointers
    for(size_t i = 0; i < bri->record_count; ++i) {
        bri->records[i].read_name.ptr = bri->readnames + bri->records[i].read_name.offset;
#ifdef BRI_INDEX_DEBUG
        REprintf("[bri-load] record %zu %s %zu\n", i, bri->records[i].read_name.ptr, bri->records[i].file_offset);
#endif
    }

    fclose(fp);
    free(index_fn);
    return bri;
}

//
// Getopt
//
// enum {
//     OPT_HELP = 1,
// };
//
// static const char* shortopts = ":i:v"; // placeholder
// static const struct option longopts[] = {
//     { "help",                      no_argument,       NULL, OPT_HELP },
//     { "index",               required_argument,       NULL,      'i' },
//     { "verbose",                   no_argument,       NULL,      'v' },
//     { NULL, 0, NULL, 0 }
// };

//
void print_usage_index()
{
    Rf_error("usage: bri index [-v] [-i <index_filename.bri>] <input.bam>\n");
}

//
// int bam_read_idx_index_main(int argc, char** argv)
// {
//     char* output_bri = NULL;
//
//     int die = 0;
//     for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
//         switch (c) {
//             case OPT_HELP:
//                 print_usage_index();
//             case 'i':
//                 output_bri = optarg;
//                 break;
//             case 'v':
//                 verbose = 1;
//                 break;
//         }
//     }
//
//     if (argc - optind < 1) {
//         Rf_error("bri index: not enough arguments\n");
//         die = 1;
//     }
//
//     if(die) {
//         print_usage_index();
//     }
//
//     char* input_bam = argv[optind++];
//     bam_read_idx_build(input_bam, output_bri);
//
//     return 0;
// }
