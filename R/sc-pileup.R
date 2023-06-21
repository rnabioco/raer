#' Pileup sites per cell
#'
#' @description This function performs a pileup operation at specified
#'   sites, returning counts for Reference (e.g. unedited) or Alternate (e.g.
#'   edited) bases. `pileup_cells` can process either a 10x Genomic's style
#'   library, from a aligned BAM file containing a tag with a cell-barcode and a
#'   tag with a UMI, or smart-seq style libraries without cell-barcodes.
#'
#'   The `sites` parameter specifies sites to pileup. This must be a [GRanges]
#'   object with 1 base intervals, a strand (+ or -), and supplemented with
#'   metadata columns named `REF` and `ALT` containing the reference and
#'   alternate base to query. See examples for an example format.
#'
#'   At each site, bases from overlapping reads will be examined, and counts of
#'   each ref and alt base enumerated for each cell-barcode present. A single
#'   base will be counted once for each UMI sequence present in each cell.
#'
#' @param bamfiles a path to a BAM file (for 10x libraries), or a vector of paths
#' to BAM files (smart-seq2). Can be supplied as a character vector, [BamFile], or
#' [BamFileList].
#' @param sites a GRanges object containing sites to process. See examples for
#'   valid formatting.
#' @param output_directory Output directory for output matrix files. The directory
#' will be generated if it doesn't exist.
#' @param chroms A character vector of chromosomes to process. If supplied, only
#'   sites present in the listed chromosomes will be processed
#' @param cell_barcodes A character vector of single cell barcodes to process. If
#' processing multiple BAM files (e.g. smart-seq-2), provide a character vector
#' of unique identifiers for each input BAM, to name each BAM file in the output files.
#' @param param object of class [FilterParam()] which specify various filters to
#'   apply to reads and sites during pileup. Note that the `min_variant_reads`
#'   parameter, if set > 0, specifies the number of variant reads at a site
#'   required in order to report a site. E.g. if set to 2, then at least 2 reads
#'   (from any cell) must have a variant in order to report the site. The
#'   default of 0 reports all sites present in the `sites` object.
#' @param umi_tag tag in BAM containing the UMI sequence
#' @param cb_tag tag in BAM containing the cell-barcode sequence
#' @param paired_end set to TRUE if data is paired-end to prevent double counting
#' overlapping read pairs.
#' @param return_sce if `TRUE`, data is returned as a SingleCellExperiment, if
#'   `FALSE` a character vector of the output files, specified by
#'   `outfile_prefix`, will be returned.
#' @param verbose Display messages
#' @param BPPARAM BiocParallel instance. Parallel computation occurs across
#'   chromosomes.
#' @param ... For the generic, further arguments to pass to specific methods.
#' Unused for now.
#'
#' @returns Returns either a [SingleCellExperiment] or character vector of paths
#'   to the sparseMatrix files produced. The [SingleCellExperiment] object is
#'   populated with two assays, `nRef` and `nAlt`, which represent base counts
#'   for the reference and alternate alleles. The [rowRanges()] will contain the
#'   genomic interval for each site, along with `REF` and `ALT` columns. The
#'   rownames will be populated with the format
#'   `site_[seqnames]_[position(1-based)]_[strand]`, with `strand` being encoded
#'   as 1 = +, 2 = -, and 3 = *.
#'
#'   If `return_sce` is `FALSE` then a character vector of paths to the
#'   sparseMatrix files (`barcodes.txt.gz`, `sites.txt.gz`, `counts.mtx.gz`),
#'   will be returned. These files can be imported using [read_sparray()].
#'
#' @examples
#' library(Rsamtools)
#' library(GenomicRanges)
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#'
#' gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
#' gr$REF <- c(rep("A", 4), "T")
#' gr$ALT <- c(rep("G", 4), "C")
#'
#' cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
#' cbs <- na.omit(cbs)
#'
#' outdir <- tempdir()
#' bai <- indexBam(bam_fn)
#'
#' fp <- FilterParam(library_type = "fr-second-strand")
#' sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
#' sce
#'
#' # example of processing multiple smart-seq2 style libraries
#'
#' many_small_bams <- rep(bam_fn, 10)
#' bam_ids <- LETTERS[1:10]
#' pileup_cells(many_small_bams,
#'     sites = gr,
#'     cell_barcodes = bam_ids,
#'     cb_tag = NULL,
#'     umi_tag = NULL,
#'     paired_end = TRUE,
#'     outdir,
#'     param = fp
#' )
#'
#' unlink(c(outdir, bai))
#'
#' @importFrom GenomeInfoDb  seqinfo seqlengths
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom BiocParallel bpworkers
#' @family pileup
#' @rdname pileup_cells
#' @export
pileup_cells <- function(bamfiles,
    sites,
    cell_barcodes,
    output_directory,
    chroms = NULL,
    umi_tag = "UB",
    cb_tag = "CB",
    paired_end = FALSE,
    param = FilterParam(),
    BPPARAM = SerialParam(),
    return_sce = TRUE,
    verbose = FALSE,
    ...) {

    if(!is(bamfiles, "BamFileList")) {
        bamfiles <- BamFileList(bamfiles)
    }

    if (length(bamfiles) > 1) {
        process_nbam <- TRUE
        if (length(bamfiles) != length(cell_barcodes)) {
            msg <- paste(
                c(
                    "multiple bamfiles detected ",
                    "a character vector of equal length ",
                    "to the number of input bams must be ",
                    "supplied to cell_barcodes"
                ),
                collapse = ""
            )
            cli::cli_abort(msg)
        }
    } else {
        process_nbam <- FALSE
    }
    bf_exists <- file.exists(path(bamfiles))
    if (!all(bf_exists)) {
        missing_bams <- path(bamfiles[!bf_exists])
        cli::cli_abort("bamfile(s) not found: {missing_bams}")
    }

    missing_index <- is.na(index(bamfiles))
    if (any(missing_index)) {
        mi <- path(bamfiles[missing_index])
        cli::cli_abort("index file(s) not found for bam: {mi}")
    }

    if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE)
    }

    if (!is(sites, "GRanges")) {
        cli::cli_abort("sites provided must be a GRanges object")
    }

    ## set default bam flags if not supplied
    if (identical(param@bam_flags, Rsamtools::scanBamFlag())) {
        param@bam_flags <- Rsamtools::scanBamFlag(
            isSecondaryAlignment = FALSE,
            isSupplementaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE
        )
    }

    cell_barcodes <- cell_barcodes[!is.na(cell_barcodes)]
    if (is.null(cb_tag)) {
        cb_tag <- character()
    } else {
        check_tag(cb_tag)
    }
    if (is.null(umi_tag)) {
        umi_tag <- character()
    } else {
        check_tag(umi_tag)
    }

    contigs <- seqinfo_from_header(bamfiles[[1]])
    contig_info <- GenomeInfoDb::seqlengths(contigs)
    chroms_to_process <- names(contig_info)
    if (!is.null(chroms)) {
        missing_chroms <- setdiff(chroms, chroms_to_process)
        if (length(missing_chroms) > 0) {
            if (verbose) {
                msg <- "the following chromosomes are not present in the bamfile(s):\n{missing_chroms}"
                cli::cli_alert_warning(msg)
            }
        }
        chroms_to_process <- intersect(chroms, chroms_to_process)
    }
    chroms_to_process <- as.character(intersect(seqnames(sites), chroms_to_process))
    sites <- sites[seqnames(sites) %in% chroms_to_process, ]

    if (length(chroms_to_process) == 0) {
        cli::cli_abort("there are no shared chromosomes found in the sites provided and bam file")
    }

    fp <- .adjustParams(param, 1)
    fp <- .as.list_FilterParam(fp)

    lib_code <- fp$library_type

    event_filters <- unlist(fp[c(
        "trim_5p",
        "trim_3p",
        "splice_dist",
        "indel_dist",
        "homopolymer_len",
        "max_mismatch_type",
        "min_read_qual",
        "min_splice_overhang",
        "min_variant_reads"
    )])

    if (verbose) cli::cli_alert("Beginning pileup")
    bf <- path.expand(path(bamfiles))
    bfi <- path.expand(index(bamfiles))
    if (process_nbam) {
        # operate in parallel over each bam
        nwrkers <- BiocParallel::bpworkers(BPPARAM)
        tmp_bamfns <- chunk_vec(bf, nwrkers)
        tmp_bamidxs <- chunk_vec(bfi, nwrkers)
        tmp_cbs <- chunk_vec(cell_barcodes, nwrkers)
        ids <- seq_along(tmp_bamfns)
        tmp_plp_files <- bpmapply(get_sc_pileup,
            bamfn = tmp_bamfns,
            index = tmp_bamidxs,
            id = ids,
            barcodes = tmp_cbs,
            MoreArgs = list(
                chrom = character(),
                sites = sites,
                outfile_prefix = output_directory,
                cb_tag = cb_tag,
                umi_tag = umi_tag,
                libtype_code = lib_code,
                event_filters = event_filters,
                fp = fp,
                pe = paired_end,
                verbose = verbose
            ),
            BPPARAM = BPPARAM,
            SIMPLIFY = FALSE
        )
    } else {
        # operate in parallel over each chromosome
        tmp_plp_files <- bpmapply(get_sc_pileup,
            chrom = chroms_to_process,
            id = chroms_to_process,
            MoreArgs = list(
                bamfn = bf,
                index = bfi,
                sites = sites,
                barcodes = cell_barcodes,
                outfile_prefix = output_directory,
                cb_tag = cb_tag,
                umi_tag = umi_tag,
                libtype_code = lib_code,
                event_filters = event_filters,
                fp = fp,
                pe = paired_end,
                verbose = verbose
            ),
            BPPARAM = BPPARAM,
            SIMPLIFY = FALSE
        )
    }

    bpstop(BPPARAM)

    if (verbose) cli::cli_alert("pileup completed, binding matrices")

    on.exit(unlink(unlist(tmp_plp_files)))
    tmp_plp_files <- Filter(function(x) length(x) > 0, tmp_plp_files)

    sces <- lapply(tmp_plp_files, function(x) {
        read_sparray(x[1], x[2], x[3], site_format = "index")
    })

    sces <- Filter(function(x) length(x) > 0, sces)

    outfns <- c("counts.mtx.gz", "sites.txt.gz", "barcodes.txt.gz")
    outfns <- file.path(output_directory, outfns)
    names(outfns) <- c("counts", "sites", "bc")

    if (length(sces) == 0) {
        warning("no sites reported in output")
        if (return_sce) {
            res <- SingleCellExperiment::SingleCellExperiment()
        } else {
            res <- outfns
        }
        return(res)
    }

    sce <- do.call(rbind, sces)
    ridx <- as.integer(rownames(sce))
    rowRanges(sce) <- sites[ridx]
    rownames(sce) <- site_names(rowRanges(sce))

    write_sparray(sce, outfns["counts"], outfns["sites"], outfns["bc"], ridx)

    if (return_sce) {
        res <- sce
    } else {
        res <- outfns
    }
    res
}

#' Read sparseMatrix produced by pileup_cells()
#'
#' @description Read in tables produced by `pileup_cells()` which are an
#'   extension of the matrixMarket sparse matrix format to store values for more
#'   than 1 matrix.
#'
#' The .mtx.gz files are formatted with columns:
#' 1) row index (0 based)
#' 2) column index (0 based)
#' 3) values for sparseMatrix #1 (nRef)
#' 4) values for sparseMatrix #2 (nAlt)
#'
#' @param mtx_fn .mtx.gz file path
#' @param sites_fn sites.txt.gz file path
#' @param bc_fn bcs.txt.gz file path
#' @param site_format one of `coordinate` or `index`, `coordinate` will populate
#' a SingleCellExperiment with rowRanges and rownames corresponind to genomic
#' intervals, whereas `index`` will only add row indicies to the rownames.
#' @returns a `SingleCellExperiment` object populated with `nRef` and `nAlt`
#' assays.
#'
#' @examples
#' library(Rsamtools)
#' library(GenomicRanges)
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#'
#' gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
#' gr$REF <- c(rep("A", 4), "T")
#' gr$ALT <- c(rep("G", 4), "C")
#'
#' cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
#' cbs <- na.omit(cbs)
#'
#' outdir <- tempdir()
#' bai <- indexBam(bam_fn)
#'
#' fp <- FilterParam(library_type = "fr-second-strand")
#' mtx_fns <- pileup_cells(bam_fn, gr, cbs, outdir,  return_sce = FALSE)
#' sce <- read_sparray(mtx_fns[1], mtx_fns[2], mtx_fns[3])
#' sce
#'
#' unlink(c(outdir, bai))
#'
#' @importFrom data.table fread
#' @importFrom Matrix sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom R.utils gzip
#' @export
read_sparray <- function(mtx_fn, sites_fn, bc_fn,
                         site_format = c("coordinate", "index")) {

    if (!file.size(mtx_fn) > 0) {
        return(SingleCellExperiment::SingleCellExperiment())
    }

    rnames <- data.table::fread(sites_fn,
                                sep = "\t",
                                col.names = c("index", "seqnames", "start",
                                              "strand", "REF", "ALT"),
                                colClasses = c("integer", "character", "integer",
                                               "integer", "character", "character"),
                                data.table = FALSE)
    site_format <- match.arg(site_format)
    if(site_format == "index") {
        rnames <- rnames$index
    } else if (site_format == "coordinate") {
        rnames <- rnames[, -1]
        rnames$strand <- c("+", "-")[rnames$strand]
        gr <- makeGRangesFromDataFrame(rnames,
                                       end.field = "start",
                                       keep.extra.columns = TRUE)
        rnames <- site_names(gr)
    }

    cnames <- readLines(bc_fn)
    n_skip <- ifelse(startsWith(readLines(mtx_fn, n = 1), "%%%"), 3, 0)

    dt <- data.table::fread(mtx_fn,
                            sep = " ",
                            colClasses = "integer",
                            skip = n_skip,
                            header = FALSE)
    sp_mtx_names <- c("nRef", "nAlt")
    n_sp_cols <- 2 + length(sp_mtx_names)

    if (ncol(dt) != n_sp_cols) cli::cli_abort("malformed sparseMatrix")

    sp_idxs <- 3:n_sp_cols
    sps <- lapply(sp_idxs, function(x) {
        sm <- sparseMatrix(
            i = dt[[1]],
            j = dt[[2]],
            x = dt[[x]],
            dims = c(length(rnames), length(cnames)),
        )
    })
    names(sps) <- sp_mtx_names
    res <- SingleCellExperiment::SingleCellExperiment(sps)
    colnames(res) <- cnames
    if(site_format == "coordinate") {
        rowRanges(res) <- gr
    }
    rownames(res) <- rnames
    res
}

write_sparray <- function(sce, mtx_fn, sites_fn, bc_fn, r_index) {
    if(!all(c("nRef", "nAlt") %in% assayNames(sce))) {
        cli::cli_abort("missing required asssays nRef or nAlt")
    }
    nref <- assay(sce, "nRef")
    nalt <- assay(sce, "nAlt")

    if(!is(nref, 'sparseMatrix')) {
        cli::cli_abort("nRef must be a sparseMatrix")
    }

    if(!is(nalt, 'sparseMatrix')) {
        cli::cli_abort("nAlt must be a sparseMatrix")
    }

    nref_trpl <- summary(nref)
    nalt_trpl <- summary(nalt)
    conforms <- identical(dim(nref_trpl), dim(nalt_trpl))
    if(!conforms){
        cli::cli_abort("nRef and nAlt sparseMatrices triplet dimensions differ")
    }

    writeLines(
        c("%%% raer MatrixMarket-like matrix coordinate integer general",
          paste("%%% ", nref@Dim[1], nref@Dim[2], length(nref@x)),
          "%%% x y nRef nAlt"),
        gzfile(mtx_fn))

    mtx <- matrix(0L, nrow = dim(nref_trpl)[1], ncol = 4L)
    mtx <- cbind(nref_trpl, nalt = nalt_trpl$x)

    data.table::fwrite(mtx, mtx_fn,
                       append = TRUE,
                       sep = " ",
                       row.names = FALSE,
                       col.names = FALSE)

    sites <- data.frame(
        r_index,
        seqnames(sce),
        start(sce),
        as.integer(strand(sce)),
        rowData(sce)$REF,
        rowData(sce)$ALT
    )

    data.table::fwrite(sites, sites_fn,
                       sep = "\t",
                       row.names = FALSE,
                       col.names = FALSE)

    writeLines(colnames(sce), gzfile(bc_fn))
}



# Utilities -------------------------------------------------------

get_sc_pileup <- function(bamfn, index, id, sites, barcodes,
    outfile_prefix, chrom,
    umi_tag, cb_tag, libtype_code,
    event_filters, fp, pe, verbose) {
    if (length(chrom) > 0) {
        sites <- sites[seqnames(sites) %in% chrom, ]
    }

    if (length(sites) == 0) {
        return(character())
    }

    outfns <- c("counts.mtx", "sites.txt", "bcs.txt")
    plp_outfns <- file.path(outfile_prefix, paste0(id, "_", outfns))
    plp_outfns <- path.expand(plp_outfns)

    if (verbose) cli::cli_alert("working on group: {id}")
    lst <- gr_to_regions(sites)
    res <- .Call(
        ".scpileup",
        bamfn,
        index,
        chrom,
        lst,
        barcodes,
        cb_tag,
        event_filters,
        fp$min_mapq,
        fp$max_depth,
        fp$min_base_quality,
        fp$read_bqual,
        as.integer(libtype_code),
        fp$bam_flags,
        plp_outfns,
        umi_tag,
        FALSE,
        pe,
        fp$min_variant_reads
    )

    if (res < 0) cli::cli_abort("pileup failed")
    plp_outfns
}

gr_to_regions <- function(gr) {
    stopifnot(all(width(gr) == 1))

    nr <- length(gr)
    if (nr == 0) cli::cli_abort("No entries in GRanges")

    gr_nms <- c("REF", "ALT")
    if (!all(gr_nms %in% names(mcols(gr)))) {
        cli::cli_abort("GRanges must have a REF and ALT columns")
    }

    if (any(strand(gr) == "*")) {
        cli::cli_alert_warning("missing strand not found in input, coercing strand to '+'")
        strand(gr) <- "+"
    }
    gr$idx <- seq(1, nr) # is one-based index
    list(
        as.character(seqnames(gr)),
        as.integer(start(gr)),
        as.integer(strand(gr)),
        as.character(mcols(gr)$REF),
        as.character(mcols(gr)$ALT),
        as.integer(mcols(gr)$idx)
    )
}

id_to_gr <- function(x, seq_info) {
    xx <- sub("^site_", "", x)
    xx <- strsplit(xx, "_")

    xx <- lapply(xx, function(x) {
        l <- length(x)
        c(seqnames = paste(x[1:(l - 4)], collapse = ""),
             start = x[[l - 3]],
             strand = x[[l - 2]],
             ref = x[[l - 1]],
             alt = x[[l]])
    })
    xx <- do.call(rbind, xx)

    if (ncol(xx) != 5) {
        cli::cli_abort("unable to decode sites to rownames")
    }

    gr <- GRanges(
        seqnames = xx[, 1],
        IRanges(
            start = as.integer(xx[, 2]),
            width = 1L
        ),
        strand = c("+", "-")[as.integer(xx[, 3])],
        ref = xx[, 4],
        alt = xx[, 5],
        seqinfo = seq_info
    )
    names(gr) <- x
    gr
}
