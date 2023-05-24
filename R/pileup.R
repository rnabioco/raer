#' Generate base counts using pileup
#'
#' @details Multiple bam files can be processed together, with files being
#'   written for each bam file. In this mode the output regions will be
#'   consistent across all files. The min_mapq, only_keep_variants, and
#'   library_type parameters can be specified for each input files.
#'
#' @param bamfiles character vector of paths to 1 or more bam files. If named,
#' the names will be included in the colData of the [RangedSummarizedExperiment], otherwise
#' the [colData] will be populated with the basename of the bamfile.
#' @param fafile path to fasta file
#' @param sites a [GRanges] object containing regions or sites to process.
#' @param region samtools region query string (i.e. chr1:100-1000). Can be combined
#' with sites, in which case sites will be filtered to keep only sites within the
#' region.
#' @param chroms chromosomes to process, not to be used with region.
#' @param param object of class [FilterParam()] which specify various
#'   filters to apply to reads and sites during pileup.
#' @param reads if supplied a fasta file will be written with reads that pass
#'   filters and contain variants
#' @param bad_reads a textfile containing read names to exclude from pileup.
#'   Readnames should be formated as readid_1 or readid_2 or readid for paired
#'   end first read paired-end second read or single end data.
#' @param umi_tag The bam tag containing a UMI sequence. If supplied, multiple
#'   reads with the same UMI sequence will only be counted once per position.
#' @param return_data if `TRUE`, data is returned as a RangedSummarizedExperiment,
#'   if `FALSE` a character vector of tabix-index files, specified by
#'   `outfile_prefix`, will be returned.
#' @param outfile_prefix Output prefix for tabix indexed files. If `NULL`, no
#'   files will be produced.
#' @param BPPARAM A [BiocParallel] class to control parallel execution. Parallel
#'   processing occurs per chromosome, so is disabled when run on a single
#'   region.
#' @param verbose if TRUE, then report progress and warnings.
#'
#' @returns A [RangedSummarizedExperiment] or a
#'   vector of the output tabix'ed file names if `return_data` is FALSE.
#'   The [RangedSummarizedExperiment] object is populated with multiple assays:
#'   * `ALT`:  Alternate base(s) found at each position
#'   * `nRef`: # of reads supporting the reference base
#'   * `nAlt`: # of reads supporting an alternate base
#'   * `nA`: # of reads with A
#'   * `nT`: # of reads with T
#'   * `nC`: # of reads with C
#'   * `nG`: # of reads with G
#'
#'   The [rowRanges()] contains the genomic interval for each site, along with:
#'   * `REF`: The reference base
#'   * `rpbz`: Mann-Whitney U test of Read Position Bias from bcftools,
#'     extreme negative or positive values indicate more bias.
#'   * `vdb`: Variant Distance Bias for filtering splice-site artefacts from
#'     bcftools, lower values indicate more bias.
#'   * `sor` Strand Odds Ratio Score, strand bias estimated by the Symmetric
#'     Odds Ratio test, based on GATK code. Higher values indicate more bias.
#'
#'   The rownames will be populated with the format
#'   `site_[seqnames]_[position(1-based)]_[strand]`, with `strand` being encoded
#'   as 1 = +, 2 = -, and 3 = *.
#'
#'   If `return_data` is `FALSE` then a character vector of paths to the
#'   pileup files will be returned. These files can be imported using
#'   `read_pileup()`.
#'
#' @examples
#' library(SummarizedExperiment)
#' bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
#' bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
#' fafn <- raer_example("human.fasta")
#'
#' rse <- pileup_sites(bamfn, fafn)
#'
#' fp <- FilterParam(only_keep_variants = TRUE, min_depth = 55)
#' pileup_sites(bamfn, fafn, param = fp)
#'
#'
#' # using multiple bam files
#'
#' bams <- rep(c(bamfn, bam2fn), each = 3)
#' sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
#' names(bams) <- sample_ids
#'
#' fp <- FilterParam(only_keep_variants = TRUE)
#' rse <- pileup_sites(bams, fafn, param = fp)
#' rse
#'
#' rse$condition <- substr(rse$sample, 1, 2)
#' assays(rse)
#'
#' colData(rse)
#'
#' rowRanges(rse)
#'
#' # specifying regions to query use GRanges object
#' sites <- rowRanges(rse)
#' rse <- pileup_sites(bams, fafn, sites = sites)
#' rse
#'
#' rse <- pileup_sites(bams, fafn, chroms = c("SPCS3", "DHFR"))
#' rse
#'
#' rse <- pileup_sites(bams, fafn, region = "DHFR:100-101")
#' rse
#' @importFrom Rsamtools bgzip indexTabix TabixFile scanTabix scanFaIndex
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels seqinfo seqlengths
#' @importFrom BiocParallel SerialParam bpstop bplapply
#'
#' @family pileup
#'
#' @export
pileup_sites <- function(bamfiles,
    fafile,
    sites = NULL,
    region = NULL,
    chroms = NULL,
    param = FilterParam(),
    outfile_prefix = NULL,
    reads = NULL,
    return_data = TRUE,
    BPPARAM = SerialParam(),
    bad_reads = NULL,
    umi_tag = NULL,
    verbose = FALSE) {
    if (is.null(names(bamfiles))) {
        sample_ids <- basename(bamfiles)
    } else {
        sample_ids <- names(bamfiles)
    }

    bamfiles <- path.expand(bamfiles)
    fafile <- path.expand(fafile)
    n_files <- length(bamfiles)

    if (!is.null(sites)) {
        if (!is(sites, "GRanges")) {
            cl <- class(sites)
            cli::cli_abort("invalid object passed to sited, expecting GRanges found {cl}")
        }
    }

    if (!all(file.exists(bamfiles))) {
        missing_bams <- bamfiles[!file.exists(bamfiles)]
        cli::cli_abort("bamfile(s) not found: {missing_bams}")
    }

    if (!file.exists(fafile)) {
        cli::cli_abort("fasta file not found: {fafile}")
    }

    if (is.null(outfile_prefix)) {
        if (!return_data) {
            cli::cli_abort("an outfile_prefix must be supplied if data is written to files")
        }
        in_memory <- TRUE
        outfiles <- character()
    } else {
        in_memory <- FALSE
        outfiles <- c(
            paste0(outfile_prefix, ".sites.txt"),
            paste0(outfile_prefix, "_", seq_len(n_files), ".plp")
        )
        if (!dir.exists(dirname(outfile_prefix))) {
            dir.create(dirname(outfile_prefix), recursive = TRUE)
        }
        outfiles <- path.expand(outfiles)
        if (length(outfiles) != n_files + 1) {
            cli::cli_abort("# of outfiles does not match # of bam input files: {outfiles}")
        }
        # remove files if exist, to avoid appending to existing files
        unlink(outfiles)

        rdat_outfn <- outfiles[1]
        plp_outfns <- outfiles[2:length(outfiles)]
    }

    contigs <- GenomeInfoDb::seqinfo(Rsamtools::BamFile(bamfiles[1]))
    contig_info <- GenomeInfoDb::seqlengths(contigs)

    if (is.null(chroms)) {
        chroms_to_process <- names(contig_info)
    } else {
        chroms_to_process <- chroms
    }

    if (is.null(region)) {
        if (!is.null(chroms)) {
            if (length(chroms) == 1) {
                region <- chroms
            }
        }
    }

    missing_chroms <- chroms_to_process[!chroms_to_process %in% names(contig_info)]

    if (length(missing_chroms) > 0) {
        if (verbose) {
            msg <- paste(missing_chroms, collapse = "\n")
            cli::cli_warn(c(
                "the following chromosomes are not present in the bamfile(s):",
                "{msg}"
            ))
        }
        chroms_to_process <- setdiff(chroms_to_process, missing_chroms)
    }

    chroms_in_fa <- seqnames(Rsamtools::scanFaIndex(fafile))
    missing_chroms <- chroms_to_process[!chroms_to_process %in% levels(chroms_in_fa)]

    if (length(missing_chroms) > 0) {
        if (verbose) {
            msg <- paste(missing_chroms, collapse = "\n")
            cli::cli_warn(c(
                "the following chromosomes are not present in the fasta file:",
                "msg"
            ))
        }
        chroms_to_process <- setdiff(chroms_to_process, missing_chroms)
    }

    chroms_to_process <-
        chroms_to_process[order(match(chroms_to_process, names(contig_info)))]

    if (!is.null(sites)) {
        sites <- sites[seqnames(sites) %in% chroms_to_process]
        sites <- gr_to_cregions(sites)
    }

    if (length(chroms_to_process) == 0) {
        cli::cli_abort("No chromosomes requested are found in bam file")
    }

    if (!is.null(reads)) {
        if (!is.character(reads) | length(reads) != 1) {
            cli::cli_abort("reads must be a character vector of length 1")
        }
    } else {
        reads <- character()
    }

    if (!is.null(bad_reads)) {
        if (!is.character(bad_reads) | length(bad_reads) != 1) {
            cli::cli_abort("bad_reads must be a character vector of length 1")
        }
    } else {
        bad_reads <- character()
    }

    if (!is.null(umi_tag)) {
        if (nchar(umi_tag) != 2) {
            cli::cli_abort("umi_tag must be a character(1) with nchar of 2 ")
        }
    } else {
        umi_tag <- character()
    }

    ## set default bam flags if not supplied
    if (identical(param@bam_flags, Rsamtools::scanBamFlag())) {
        param@bam_flags <- Rsamtools::scanBamFlag(
            isSecondaryAlignment = FALSE,
            isNotPassingQualityCont = FALSE,
            isDuplicate = FALSE,
            isSupplementaryAlignment = FALSE
        )
    }

    fp <- .adjustParams(param, n_files)
    fp <- .c_args_FilterParam(fp)

    if (!is.null(region)) {
        chroms_to_process <- region
    } else {
        region <- character()
    }

    if (is(BPPARAM, "SerialParam") || length(chroms_to_process) == 1) {
        start_time <- Sys.time()
        if (length(chroms_to_process) > 1 && is.null(sites)) {
            to_process <- contig_info[chroms_to_process]
            sites <- GRanges(paste0(names(to_process), ":", 1, "-", to_process))
            sites <- gr_to_cregions(sites)
        }
        res <- .Call(
            ".pileup", bamfiles, as.integer(n_files), fafile, region,
            sites, fp[["int_args"]], fp[["numeric_args"]], fp[["lgl_args"]],
            fp[["library_type"]], fp[["only_keep_variants"]], fp[["min_mapq"]],
            in_memory, outfiles, reads, bad_reads, umi_tag
        )

        if (!in_memory) {
            if (res != 0) {
                cli::cli_abort("Error occured during pileup")
            }
        } else {
            rdat <- res[[1]]
            res <- res[2:length(res)]
            res <- lists_to_grs(res, contigs)
            res <- merge_pileups(res,
                sample_names = sample_ids
            )
            rowData(res) <- cbind(rowData(res), rdat)
        }

        if (verbose) {
            time_elapsed <- prettyNum(Sys.time() - start_time, digits = 3)
            cli::cli_alert_success("Completed pileup in {time_elapsed}")
        }
    } else {
        res <- bplapply(chroms_to_process, function(ctig) {
            start_time <- Sys.time()
            if (length(outfiles) > 0) {
                tmp_outfiles <- unlist(lapply(
                    seq_along(outfiles),
                    function(x) tempfile()
                ))
                tmp_rdatfile <- tmp_outfiles[1]
                tmp_plpfiles <- tmp_outfiles[2:length(tmp_outfiles)]
                fn_lst <- list(
                    contig = ctig,
                    bam_fn = bamfiles,
                    tmp_plpfns = tmp_plpfiles,
                    tmp_rdatfn = tmp_rdatfile
                )
            } else {
                tmp_outfiles <- character()
            }
            if (is.null(ctig)) ctig <- character()
            res <- .Call(
                ".pileup", bamfiles, as.integer(n_files), fafile, ctig,
                sites, fp[["int_args"]], fp[["numeric_args"]],
                fp[["lgl_args"]], fp[["library_type"]], fp[["only_keep_variants"]],
                fp[["min_mapq"]], in_memory, tmp_outfiles, reads, bad_reads, umi_tag
            )

            if (!in_memory) {
                if (res != 0) {
                    cli::cli_abort("Error occured during pileup")
                }
                res <- fn_lst
            } else {
                rdat <- res[[1]]
                res <- res[2:length(res)]
                res <- lists_to_grs(res, contigs)
                res <- list(rdat = rdat, plps = res)
            }

            if (verbose) {
                time_elapsed <- prettyNum(Sys.time() - start_time, digits = 3)
                cli::cli_alert_success("Completed pileup on {ctig} in {time_elapsed}")
            }
            res
        }, BPPARAM = BPPARAM)
        bpstop(BPPARAM)
        if (!in_memory) {
            for (i in seq_along(res)) {
                ra <- file.append(rdat_outfn, res[[i]]$tmp_rdatfn)
                pa <- file.append(plp_outfns, res[[i]]$tmp_plpfns)
                if (!all(c(ra, pa))) {
                    cli::cli_abort("error occured concatenating pileup files")
                }
                unlink(c(res[[i]]$tmp_rdatfn, res[[i]]$tmp_plpfns))
            }
        } else {
            res <- lapply(res, function(x) {
                se <- merge_pileups(x$plps, sample_names = sample_ids)
                rowData(se) <- cbind(rowData(se), x$rdat)
                se
            })
            res <- do.call(rbind, res)
        }
    }

    if (!in_memory) {
        if (any(file.info(outfiles)$size == 0)) {
            return(empty_plp_record())
        }

        # run_pileup writes to a (temp)file, next the file will be tabix indexed
        tbxfiles <- lapply(outfiles, function(x) {
            tbxfile <- Rsamtools::bgzip(x, overwrite = TRUE)
            idx <- Rsamtools::indexTabix(tbxfile, seq = 1, start = 2, end = 2, zeroBased = FALSE)
            tbxfile
        })

        if (!return_data) {
            unlink(outfiles)
            return(unlist(tbxfiles))
        }

        res <- lapply(tbxfiles, function(x) {
            xx <- read_pileup(x, region = NULL)
            GenomeInfoDb::seqlevels(xx) <- GenomeInfoDb::seqlevels(contigs)
            GenomeInfoDb::seqinfo(xx) <- contigs
            xx
        })

        rdat <- res[[1]]
        res <- res[2:length(res)]
        unlink(outfiles)
    }
    res
}


# IRanges/GRanges are limited to this max int
MAX_INT <- 536870912

#' Read pileup, indexed by tabix
#'
#' @description This function can read in pileup or site data tables produced by
#' pileup_sites().
#' @param tbx_fn filename
#' @param region region to read from file, samtools style
#' region specifiers are supported.
#'
#' @examples
#' bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#' plp_fn <- tempfile()
#' plp <- pileup_sites(bamfn, fafn, return_data = FALSE, outfile_prefix = plp_fn)
#'
#' # first table contains site information, second table pileup count
#' lapply(plp, function(x) head(read_pileup(x)))
#' unlink(c(plp, plp_fn, paste0(plp, ".tbi")))
#'
#' @return `GenomicRanges::GRanges` containing positions of editing sites and pileup
#' data.
#' @importFrom data.table fread
#' @export
read_pileup <- function(tbx_fn, region = NULL) {
    tbx <- Rsamtools::TabixFile(tbx_fn)

    # using Rsamtools read in tabix file
    # note that file is read in as a list of character vectors
    # consider using our own read_tabix function if this is a bottleneck
    if (!is.null(region)) {
        ivl_vals <- get_region(region)
        # note that samtools will return a larger INT than IRANGES can handle
        # if no ranges are supplied in the region
        ivl_end <- min(MAX_INT, ivl_vals$end)
        params <- GenomicRanges::GRanges(
            ivl_vals$chrom,
            IRanges::IRanges(
                start = ivl_vals$start + 1,
                end = ivl_end
            )
        )
        tbx_vals <- Rsamtools::scanTabix(tbx, param = params)[[1]]
    } else {
        tbx_vals <- Rsamtools::scanTabix(tbx)[[1]]
    }

    # quick method to convert vector of character strings into data.frame
    if (length(tbx_vals) == 1) {
        nc <- lengths(regmatches(tbx_vals, gregexpr("\t", tbx_vals, fixed = TRUE)))

        # handle length 1 character vectors, which will not work with fread
        from <- data.frame(t(strsplit(tbx_vals, "\t")[[1]]))
        colnames(from) <- paste0("V", seq_len(ncol(from)))
        if(nc == 7) {
            from[c(2, 5:nc)] <- as.numeric(from[c(2, 5:nc)])
        } else {
            from[c(2, 6:nc)] <- as.numeric(from[c(2, 6:nc)])
        }
    } else {
        from <- data.table::fread(
            text = tbx_vals,
            stringsAsFactors = FALSE,
            data.table = FALSE,
            showProgress = FALSE,
            sep = "\t"
        )
    }

    plp_cols <- c("ALT", "nRef", "nAlt", "nA", "nT", "nC", "nG", "nN", "nX")
    rowData_cols <- c("rpbz", "vdb", "sor")
    if (ncol(from) ==  (4 + length(rowData_cols))) {
        cols <- rowData_cols
    } else if (ncol(from) == (4 + length(plp_cols))) {
        cols <- plp_cols
    } else {
        cli::cli_abort("unknown pileup file format")
    }

    colnames(from)[4:ncol(from)] <- c("REF", cols)

    GenomicRanges::GRanges(
        seqnames = from$V1,
        ranges = IRanges::IRanges(
            start = from$V2,
            end = from$V2
        ),
        strand = from$V3,
        from[, 4:ncol(from)]
    )
}

# parses a region string (e.g. "chr1:456-567" or "chr1")
# using htslib into chrom, start, and end (0-based start)
# returns a list with chrom, start, and end entries
# if only a chromosome is supplied then start is set
# to 0 and end is set to 2^31 - 1
get_region <- function(region) {
    if (!is.character(region) || length(region) != 1) {
        stop("region must be character(1)")
    }
    .Call(".get_region", region)
}

# generate empty pileup record
# idea from @user2462304 https://stackoverflow.com/a/48180979/6276041
empty_plp_record <- function() {
    col_types <- list(
        REF = character(),
        ALT = character(),
        nRef = integer(),
        nAlt = integer(),
        nA = integer(),
        nT = integer(),
        nC = integer(),
        nG = integer(),
        nN = integer(),
        nX = integer()
    )
    df <- do.call(data.frame, col_types)
    gr <- GRanges(c(seqnames = NULL, ranges = NULL, strand = NULL))
    mcols(gr) <- df
    gr
}


.adjust_arg_length <- function(obj, name, len) {
    if (length(slot(obj, name)) != len) {
        if (length(slot(obj, name)) == 1) {
            slot(obj, name) <- rep(slot(obj, name), len)
        } else {
            stop(
                "%s requires either 1 value, or individual values,",
                "for all input bamfiles", slot
            )
        }
    }
    slot(obj, name)
}
## Check validity and adjust
.adjustParams <- function(filterParam, nFiles) {
    if (!inherits(filterParam, "FilterParam")) {
        stop(
            "'filterParam' must inherit from 'FilterParam', got '%s'",
            class(filterParam)
        )
    }
    filterParam@min_mapq <- .adjust_arg_length(
        filterParam,
        "min_mapq",
        nFiles
    )
    filterParam@only_keep_variants <- .adjust_arg_length(
        filterParam,
        "only_keep_variants",
        nFiles
    )
    filterParam@library_type <- .adjust_arg_length(
        filterParam,
        "library_type",
        nFiles
    )
    filterParam
}


#' @importFrom methods slot slot<- slotNames
.FilterParam <- setClass(
    "FilterParam",
    slots = c(
        max_depth = "integer",
        min_depth = "integer",
        min_base_quality = "integer",
        min_mapq = "integer", # variable length
        library_type = "integer", # variable length
        bam_flags = "integer", # length 2
        trim_5p = "integer",
        trim_3p = "integer",
        indel_dist = "integer",
        splice_dist = "integer",
        min_splice_overhang = "integer",
        homopolymer_len = "integer",
        max_mismatch_type = "integer", # length 2
        min_variant_reads = "integer",
        only_keep_variants = "logical", # variable length
        report_multiallelic = "logical",
        ftrim_5p = "numeric",
        ftrim_3p = "numeric",
        read_bqual = "numeric", # length 2
        min_allelic_freq = "numeric"
    )
)

setMethod(show, "FilterParam", function(object) {
    cat("class: ", class(object), "\n")
    values <- lapply(slotNames(object), slot, object = object)
    info <- paste(slotNames(object), values, sep = ": ", collapse = "; ")
    cat(strwrap(info, exdent = 2), sep = "\n")
})


.encode_libtype <- function(library_type = c("genomic-unstranded",
                                             "fr-first-strand",
                                             "fr-second-strand",
                                             "unstranded"),
                            n_files) {
    # encode libtype as integer
    # 0 = genomic-unstranded  all reads on + strand
    # 1 = fr-first-strand     strand based on R1/antisense, R2/sense
    # 2 = fr-second-strand    strand based on R1/sense, R2/antisense
    # 3 = unstranded          strand based on alignment
    lib_values <- c(
        "genomic-unstranded",
        "fr-first-strand",
        "fr-second-strand",
        "unstranded"
    )
    lib_code <- match(library_type, lib_values)
    if (any(is.na(lib_code))) {
        stop("library_type must be one of :", paste(lib_values, collapse = " "))
    } else {
        lib_code <- lib_code - 1
    }

    as.integer(lib_code)
}

.c_args_FilterParam <- function(x, ...) {
    slotnames <- slotNames(x)
    names(slotnames) <- slotnames
    fp <- lapply(slotnames, slot, object = x)

    # consistent length args are populated into vectors
    # note that unlisting will increase vector size greather than number of args
    # (e.g. bam_flags will add 2 to vector)
    int_args <- unlist(fp[c(
        "max_depth",
        "min_depth",
        "min_base_quality",
        "trim_5p",
        "trim_3p",
        "indel_dist",
        "splice_dist",
        "min_splice_overhang",
        "homopolymer_len",
        "min_variant_reads",
        "max_mismatch_type", # length 2
        "bam_flags" # length 2
    )])

    numeric_args <- unlist(fp[c(
        "ftrim_5p",
        "ftrim_3p",
        "min_allelic_freq",
        "read_bqual" # length 2
    )])

    lgl_args <- unlist(fp[c(
        "report_multiallelic"
    )])
    # variable length args
    # passed as separate args to c fxns

    list(
        int_args = int_args,
        numeric_args = numeric_args,
        lgl_args = lgl_args,
        library_type = fp[["library_type"]],
        only_keep_variants = fp[["only_keep_variants"]],
        min_mapq = fp[["min_mapq"]]
    )
}

.as.list_FilterParam <- function(x, ...) {
    slotnames <- slotNames(x)
    names(slotnames) <- slotnames
    lapply(slotnames, slot, object = x)
}

#' @param min_depth min read depth needed to report site
#' @param max_depth maximum read depth considered at each site
#' @param min_base_quality min base quality score to consider read for pileup
#' @param min_mapq minimum required MAPQ score, can be a vector of values
#' for each bam file
#' @param library_type read orientation, one of `fr-first-strand`,
#' `fr-second-strand`, `unstranded`, and `genomic-unstranded`. Can supply as a vector to specify for each
#' input bam. Unstranded library type will be reported based on read alignment.
#' genomic-unstranded will report all variants w.r.t the + strand.
#' @param only_keep_variants if TRUE, then only variant sites will be reported
#' (FALSE by default), can be a vector for each input bamfile
#' @param bam_flags bam flags to filter or keep, use [Rsamtools::scanBamFlag()]
#'   to generate.
#' @param trim_5p Bases to trim from 5' end of read alignments
#' @param trim_3p Bases to trim from 3' end of read alignments
#' @param ftrim_5p Fraction of bases to trim from 5' end of read alignments
#' @param ftrim_3p Fraction of bases to trim from 3' end of read alignments
#' @param splice_dist Exclude read if site occurs within given
#' distance from splicing event in the read
#' @param min_splice_overhang Exclude read if site is located adjacent to splice
#' site with an overhang of less than given length.
#' @param indel_dist Exclude read if site occurs within given
#' distance from indel event in the read
#' @param homopolymer_len Exclude site if occurs within homopolymer of given
#' length
#' @param max_mismatch_type Exclude read if it has X different mismatch types
#' (e.g A-to-G, G-to-C, C-to-G, is 3 mismatch types) or Y # of mismatches,
#' must be supplied as a integer vector of length 2. e.g.
#' c(X, Y).
#' @param read_bqual Exclude read if more than X percent of the bases have
#' base qualities less than Y. Numeric vector of length 2. e.g. c(0.25, 20)
#' @param min_variant_reads Required number of reads containing a variant for a site
#' to be reported. Calculated per bam file, such that if 1 bam file has >= min_variant_reads,
#' then the site will be reported.
#' @param min_allelic_freq minimum allelic frequency required for a variant to be
#' reported in ALT assays.
#' @param report_multiallelic if TRUE, report sites with multiple variants passing
#' filters. If FALSE, site will not be reported.
#'
#' @rdname pileup_sites
#' @export
FilterParam <-
    function(max_depth = 1e4, min_depth = 1L, min_base_quality = 20L,
    min_mapq = 0L, library_type = "fr-first-strand",
    bam_flags = NULL, only_keep_variants = FALSE,
    trim_5p = 0L, trim_3p = 0L, ftrim_5p = 0, ftrim_3p = 0,
    indel_dist = 0L, splice_dist = 0L, min_splice_overhang = 0L,
    homopolymer_len = 0L,
    max_mismatch_type = c(0L, 0L), read_bqual = c(0.0, 0.0),
    min_variant_reads = 0L, min_allelic_freq = 0,
    report_multiallelic = TRUE) {
        stopifnot(isSingleNumber(max_depth))
        stopifnot(isSingleNumber(min_base_quality))
        stopifnot(isSingleNumber(min_depth))
        stopifnot(isSingleNumber(trim_5p))
        stopifnot(isSingleNumber(trim_3p))
        stopifnot(isSingleNumber(indel_dist))
        stopifnot(isSingleNumber(splice_dist))
        stopifnot(isSingleNumber(homopolymer_len))
        stopifnot(isSingleNumber(min_splice_overhang))
        stopifnot(isSingleNumber(min_variant_reads))
        stopifnot(isSingleNumber(ftrim_5p))
        stopifnot(isSingleNumber(ftrim_3p))
        stopifnot(isSingleNumber(min_allelic_freq))

        max_depth <- as.integer(max_depth)
        min_base_quality <- as.integer(min_base_quality)

        min_depth <- as.integer(min_depth)
        trim_5p <- as.integer(trim_5p)
        trim_3p <- as.integer(trim_3p)
        indel_dist <- as.integer(indel_dist)
        splice_dist <- as.integer(splice_dist)
        min_mapq <- as.integer(min_mapq)
        homopolymer_len <- as.integer(homopolymer_len)
        max_mismatch_type <- as.integer(max_mismatch_type)
        read_bqual <- as.numeric(read_bqual)
        min_splice_overhang <- as.integer(min_splice_overhang)
        min_variant_reads <- as.integer(min_variant_reads)
        ftrim_5p <- as.numeric(ftrim_5p)
        ftrim_3p <- as.numeric(ftrim_3p)
        min_allelic_freq <- as.numeric(min_allelic_freq)

        stopifnot(ftrim_5p >= 0 && ftrim_5p <= 1)
        stopifnot(ftrim_3p >= 0 && ftrim_3p <= 1)

        stopifnot(length(max_mismatch_type) == 2 && !any(is.na(max_mismatch_type)))
        stopifnot(length(read_bqual) == 2 && !any(is.na(read_bqual)))
        stopifnot(isTRUEorFALSE(report_multiallelic))

        # variable length depending on n_files
        stopifnot(is.character(library_type))
        stopifnot(is.integer(min_mapq))
        stopifnot(is.logical(only_keep_variants))

        if (is.null(bam_flags)) {
            # defaults to allowing all reads
            bam_flags <- Rsamtools::scanBamFlag()
        } else {
            if (length(bam_flags) != 2 || !all(names(bam_flags) == c("keep0", "keep1"))) {
                stop("bam_flags must be generated using Rsamtools::scanBamFlag()")
            }
        }

        library_type <- .encode_libtype(library_type)

        ## creation
        .FilterParam(
            max_depth = max_depth, min_base_quality = min_base_quality,
            min_mapq = min_mapq, min_depth = min_depth,
            library_type = library_type, only_keep_variants = only_keep_variants,
            trim_5p = trim_5p, trim_3p = trim_3p, indel_dist = indel_dist,
            splice_dist = splice_dist, homopolymer_len = homopolymer_len,
            max_mismatch_type = max_mismatch_type, read_bqual = read_bqual,
            min_splice_overhang = min_splice_overhang,
            min_variant_reads = min_variant_reads,
            ftrim_5p = ftrim_5p, ftrim_3p = ftrim_3p,
            min_allelic_freq = min_allelic_freq, report_multiallelic = report_multiallelic,
            bam_flags = bam_flags
        )
    }

PILEUP_COLS <- c(
    "seqnames",
    "pos",
    "strand",
    "REF",
    "ALT",
    "nRef",
    "nAlt",
    "nA",
    "nT",
    "nC",
    "nG",
    "nN",
    "nX"
)

# Utilities -------------------------------------------------


# convert list of lists to list of grs
lists_to_grs <- function(x, seqinfo = NULL) {
    mc_cols <- setdiff(PILEUP_COLS, c("seqnames", "pos", "strand"))
    lapply(x, function(mc) {
        GRanges(
            seqnames = mc$seqname,
            ranges = IRanges(
                start = mc$pos,
                width = 1L
            ),
            strand = mc$strand,
            mc[mc_cols],
            seqinfo = seqinfo
        )
    })
}

gr_to_cregions <- function(gr) {
    nr <- length(gr)
    if (nr == 0) cli::cli_abort("No entries in GRanges")
    list(
        as.character(seqnames(gr)),
        as.integer(start(gr)),
        as.integer(end(gr))
    )
}


#' Merge pileups into a RangedSummarizedExperiment
#'
#' @description Create a `RangedSummarizedExperiment` from a single or list of
#'   pileups (e.g., for different samples) generated by `pileup_sites()`.
#'
#' @param plps results from running [pileup_sites()], can be one result, a list of
#'   results, or a named list of results. If a named list is given, the colData
#'   will be named using the names in the list.
#' @param assay_cols character vector of columns to store as assays
#' @param sample_names A list of names to be added to the SE object. If no
#'   sample names are given and plps is not a named list, then default names (ie
#'   sample_1, sample_2, ..., sample_n) will be given and a warning will be
#'   printed.
#' @param fill_na Numeric value to replace NAs in numeric matrices Should only
#'   be used when plps were computed independently with a min_nucleotide_count =
#'   1, otherwise sites may be set to 0, although they may have coverage > 0 but
#'   less than the min_nucleotide_count parameter.
#' @param verbose print information on progress
#'
#' @return `RangedSummarizedExperiment` populated with assays specified in
#'   `assay_cols`.
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom IRanges extractList
#' @keywords internal
#' @noRd
merge_pileups <- function(plps,
    assay_cols = c("ALT", "nRef", "nAlt", "nA", "nT", "nC", "nG"),
    sample_names = NULL,
    fill_na = NULL,
    verbose = FALSE) {
    if (!is.list(plps)) {
        plps <- list(plps)
    }
    # Checks for sample names
    if (is.null(sample_names)) {
        if (is.null(names(plps))) {
            sample_names <- paste0("sample_", seq_along(plps))
        } else {
            sample_names <- names(plps)
        }
    } else {
        if (length(plps) != length(sample_names)) {
            cli::cli_abort(c(
                "You must provide the same number of sample names as pileup results",
                "You supplied {length(plps)} pileup results but supplied ",
                "{length(sample_names} sample names."
            ))
        }
        names(plps) <- sample_names
    }

    stopifnot(all(assay_cols %in% colnames(mcols(plps[[1]]))))

    # determine all ranges across plps,
    # store Ref nt in mcols (should be invariant)
    all_ranges <- unlist(as(plps, "GRangesList"), use.names = FALSE)
    ref_nt <- mcols(all_ranges)[["REF"]]
    mcols(all_ranges) <- NULL
    mcols(all_ranges)$REF <- ref_nt
    all_ranges <- sort(unique(all_ranges), ignore.strand = TRUE)
    assay_cols <- setdiff(assay_cols, "REF")

    plp_assays <- fill_matrices(plps, assay_cols, all_ranges, verbose)

    se <- SummarizedExperiment(
        assays = plp_assays,
        rowRanges = all_ranges,
        colData = DataFrame(sample = sample_names)
    )
    colnames(se) <- se$sample
    rownames(se) <- site_names(rowRanges(se))

    if (!is.null(fill_na)) {
        assays(se) <- lapply(assays(se), function(x) {
            if (is.numeric(x)) x[is.na(x)] <- fill_na
            x
        })
    }
    se
}


# standardize and bind matrices from pileups
# into dense matrices
fill_matrices <- function(plps, assays, gr, verbose = TRUE) {
    if (length(assays) == 0) {
        return(NULL)
    }

    plps <- as(plps, "GRangesList")
    plp_assays <- lapply(assays, function(x) {
        matrix(NA,
            nrow = NROW(gr),
            ncol = length(plps)
        )
    })
    names(plp_assays) <- assays

    for (i in seq_along(plps)) {
        if (verbose) {
            if (i %% 128 == 0) {
                message("processed ", i, " pileups")
            }
        }
        hits <- findOverlaps(gr, plps[[i]], type = "equal")
        hits <- queryHits(hits)
        for (j in seq_along(plp_assays)) {
            id <- names(plp_assays)[j]
            plp_assays[[j]][hits, i] <- mcols(plps[[i]])[[id]]
        }
    }
    plp_assays
}

bind_se <- function(ses) {
    all_assays <- lapply(ses, function(x) names(assays(x)))
    assay_names <- unique(unlist(all_assays))
    stopifnot(all(unlist(lapply(
        all_assays,
        function(x) {
            length(setdiff(assay_names, x)) == 0
        }
    ))))

    assay_dat <- vector("list", length = length(assay_names))
    for (i in seq_along(assay_names)) {
        an <- assay_names[i]
        an_assays <- lapply(ses, function(x) assay(x, an))
        assay_dat[[i]] <- cbind_sparse(an_assays)
        names(assay_dat)[i] <- an
    }

    cdat <- do.call(rbind, lapply(ses, colData))

    rns <- lapply(assay_dat, rownames)
    stopifnot(all(unlist(lapply(rns[-1], identical, rns[[1]]))))

    rn <- rownames(assay_dat[[1]])
    rowrng <- unlist(GRangesList(lapply(ses, rowRanges)), use.names = FALSE)
    rowrng <- rowrng[rn, , drop = FALSE]
    SummarizedExperiment(
        assays = assay_dat,
        rowRanges = rowrng,
        colData = cdat
    )
}


# adapted from user6376297
# https://stackoverflow.com/a/56547070/6276041
# cbind sparse matrices with differing row entries
cbind_sparse <- function(mats) {
    stopifnot(all(unlist(lapply(mats, function(x) is(x, "dgCMatrix")))))
    rn <- unique(unlist(lapply(mats, rownames)))
    cn <- unique(unlist(lapply(mats, colnames)))
    nvals <- sum(unlist(lapply(mats, Matrix::nnzero)))

    if (nvals > 0) {
        x <- integer(nvals)
        i <- integer(nvals)
        j <- integer(nvals)

        init_x <- 1
        for (idx in seq_along(mats)) {
            cindnew <- match(colnames(mats[[idx]]), cn)
            rindnew <- match(rownames(mats[[idx]]), rn)
            ind <- summary(mats[[idx]])
            nval <- nrow(ind)
            end <- nval + init_x - 1
            # handle binding a zero entry matrix to a non-zero
            if (nval > 0) {
                i[init_x:end] <- rindnew[ind$i]
                j[init_x:end] <- cindnew[ind$j]
                x[init_x:end] <- ind$x
            }
            init_x <- end + 1
        }

        res <- sparseMatrix(
            i = i,
            j = j,
            x = x,
            dims = c(length(rn), length(cn)),
            dimnames = list(rn, cn)
        )
    } else {
        # handle binding a zero entry matrices
        res <- sparseMatrix(
            i = integer(0),
            p = 0,
            dims = c(length(rn), length(cn)),
            dimnames = list(rn, cn)
        )
    }
    res
}

site_names <- function(gr) {
    if (length(gr) == 0) {
        return(NULL)
    }
    paste(
        "site",
        decode(seqnames(gr)),
        decode(start(gr)),
        as.integer(strand(gr)),
        sep = "_"
    )
}
