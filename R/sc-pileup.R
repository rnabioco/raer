#' Pileup sites per cell
#'
#' @description This function will perform a pileup operation at specified
#'   sites, returning counts for Reference (e.g. unedited) or Alternate (e.g.
#'   editing) bases. Current functionality will process a 10x genomic's style
#'   library, from a aligned bam file containing a tag with a cell-barcode and a
#'   tag with a UMI.
#'
#'   The `sites` parameter specifies sites to pileup. This must be a GRanges
#'   object with 1 base intervals, a strand (+ or -), and supplemented with
#'   metadata columns named `ref` and `alt` containing the reference and
#'   alternate base to query. See examples for an example format.
#'
#'   At each site, bases from overlapped reads will be examined, and counts of
#'   each ref and alt base enumerated for each cell barcode present. A single
#'   base will be counted once for each UMI sequence present in each cell.
#'
#' @param bamfile BAM file name
#' @param sites a GRanges object containing sites to process. See examples for
#'   valid formatting.
#' @param output_directory Output directory for output files. Will be generated
#'   if it doesn't exist.
#' @param chroms A character vector of chromosomes to process, if supplied, only
#'   sites present in the listed chromosomes will be processed
#' @param cell_barcodes A character vector of single cell barcodes to process.
#' @param param object of class [FilterParam()] which specify various filters to
#'   apply to reads and sites during pileup. Note that the `min_variant_reads`
#'   parameter, if set > 0, specifies the number of variant reads at a site
#'   required in order to report a site. E.g. if set to 2, then at least 2 reads
#'   (from any cell) must have a variant in order to report the site. The
#'   default of 0 reports all sites present in the `sites` object.
#' @param umi_tag tag in bam containing the UMI sequence
#' @param cb_tag tag in bam containing the cell barcode sequence
#' @param paired_end set to TRUE if data is paired end.
#' @param return_sce if `TRUE`, data is returned as a SingleCellExperiment, if
#'   `FALSE` a character vector of the output files, specified by
#'   `outfile_prefix`, will be returned.
#' @param verbose Display messages
#' @param BPPARAM BiocParallel instance. Parallel computation occurs across
#'   chromosomes.
#'
#' @returns Returns either a `SingleCellExperiment` or character vector of paths
#'   to the files produced.
#'
#' @examples
#' library(Rsamtools)
#' library(GenomicRanges)
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#'
#' gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-","2:589:-", "2:601:-"))
#' gr$ref <- c(rep("A", 4), "T")
#' gr$alt <- c(rep("G", 4), "C")
#'
#' cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
#' cbs <- na.omit(cbs)
#'
#' outdir <- tempdir()
#' bai <- indexBam(bam_fn)
#' on.exit(unlink(outdir, bai))
#'
#' fp <- FilterParam(library_type = "fr-second-strand")
#' sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
#' sce
#'
#' @importFrom GenomeInfoDb  seqinfo seqlengths
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom BiocParallel bpworkers
#' @family pileup
#'
#' @export
pileup_cells <- function(bamfile,
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
                         verbose = FALSE) {

  if(length(bamfile) > 1) {
    process_nbam <- TRUE
    if(length(bamfile) != length(cell_barcodes)){
      msg <- paste(c("multiple bamfiles detected ",
                     "a character vector of equal length ",
                     "to the number of input bams must be ",
                     "supplied to cell_barcodes"),
                   collapse = "")
      cli::cli_abort(msg)
    }
  } else {
    process_nbam <- FALSE
  }
  if (!all(file.exists(bamfile))) {
    missing_bams <- bamfile[!file.exists(bamfile)]
    cli::cli_abort("bamfile(s) not found: {missing_bams}")
  }

  bamfile <- path.expand(bamfile)
  seq_info <- GenomeInfoDb::seqinfo(Rsamtools::BamFile(bamfile[1]))

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  if(!is(sites, "GRanges")){
    cli::cli_abort("sites provided are not formatted as GRanges")
  }

  ## set default bam flags if not supplied
  if(identical(param@bam_flags, Rsamtools::scanBamFlag())){
    param@bam_flags <- Rsamtools::scanBamFlag(
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE
    )
  }

  cell_barcodes <- cell_barcodes[!is.na(cell_barcodes)]
  if(is.null(cb_tag)){
    cb_tag = character()
  } else {
    check_tag(cb_tag)
  }
  if(is.null(umi_tag)){
    umi_tag = character()
  } else {
    check_tag(umi_tag)
  }

  contigs <- GenomeInfoDb::seqinfo(Rsamtools::BamFile(bamfile[1]))
  contig_info <- GenomeInfoDb::seqlengths(contigs)
  chroms_to_process <- names(contig_info)
  if (!is.null(chroms)) {
    missing_chroms <- setdiff(chroms, chroms_to_process)
    if (length(missing_chroms) > 0) {
      if(verbose){
        msg <- "the following chromosomes are not present in the bamfile(s):\n{missing_chroms}"
        cli::cli_alert_warning(msg)
      }
    }
    chroms_to_process <- intersect(chroms, chroms_to_process)
  }
  chroms_to_process <- as.character(intersect(seqnames(sites), chroms_to_process))
  sites <- sites[seqnames(sites) %in% chroms_to_process, ]

  if(length(chroms_to_process) == 0){
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
  if(process_nbam){
    nwrkers <- BiocParallel::bpworkers(BPPARAM)
    tmp_bamfns <- chunk_vec(bamfile, nwrkers)
    tmp_cbs <- chunk_vec(cell_barcodes, nwrkers)
    ids <- seq_along(tmp_bamfns)
    tmp_plp_files <- bpmapply(get_sc_pileup,
                              bamfn = tmp_bamfns,
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
                              SIMPLIFY = FALSE)
  } else {
    tmp_plp_files <- bpmapply(get_sc_pileup,
                              chrom = chroms_to_process,
                              id = chroms_to_process,
                              MoreArgs = list(
                                bamfn = bamfile,
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
                              SIMPLIFY = FALSE)
  }

  bpstop(BPPARAM)

  if (verbose) cli::cli_alert("pileup completed, binding matrices")

  tmp_plp_files <- Filter(function(x) length(x) > 0, tmp_plp_files)
  sps <- lapply(tmp_plp_files, function(x){
    read_sparray(x[1], x[2], x[3])
  })
  sps <- sps[!unlist(lapply(sps, is.null))]
  sp_assays <- list(nRef = lapply(sps, function(x) x[[1]]),
                    nVar = lapply(sps, function(x) x[[2]]))

  outfns <- c("barcodes.txt.gz","sites.txt.gz","counts.mtx.gz")
  outfns <- file.path(output_directory, outfns)
  names(outfns) <- c("bc", "sites", "counts")

  sp_assays <- Filter(function(x) length(x) > 0, sp_assays)
  if(length(sp_assays) == 0) {
    warning("no sites reported in output")
    if(return_sce){
      res <- SingleCellExperiment::SingleCellExperiment()
    } else {
      res <- outfns
    }
    return(res)
  }

  if(process_nbam){
    sp_assays <- lapply(sp_assays, function(x) {
      xx <- do.call(cbind, x)
      writeLines(colnames(xx), outfns["bc"])
      writeLines(rownames(xx), outfns["sites"])
      Matrix::writeMM(xx, outfns["counts"])
      xx
    })
  } else {
    sp_assays <- lapply(sp_assays, function(x) {
      xx <- do.call(rbind, x)
      writeLines(colnames(xx), outfns["bc"])
      writeLines(rownames(xx), outfns["sites"])
      Matrix::writeMM(xx, outfns["counts"])
      xx
    })
  }

  unlink(unlist(tmp_plp_files))

  if(return_sce){
    res <- SingleCellExperiment::SingleCellExperiment(sp_assays)
    res$id <- colnames(res)
    rowRanges(res) <- id_to_gr(rownames(res), seq_info)
  } else {
    res <- outfns
  }
  res
}

#' Read tables produced by pileup_cells()
#'
#' @description Read in tables produced by `pileup_cells()`, which are an
#'   extension of the matrixMarket sparse matrix format to store values for more
#'   than 1 matrix.
#'
#' The .mtx.gz files are formatted with columns:
#' 1) row index (0 based)
#' 2) column index (0 based)
#' 3) values for sparseMatrix #1 (nRef)
#' 4) values for sparseMatrix #2 (nVar)
#' N) values for sparseMatrix ... (...) ununsed for now
#'
#' @param mtx_fn .mtx.gz file path
#' @param sites_fn sites.txt.gz file path
#' @param bc_fn bcs.txt.gz file path
#'
#' @returns a list of `sparseMatrix`, of `NULL` if mtx_fn is empty
#'
#' @importFrom data.table fread
#'
#' @export
read_sparray <- function(mtx_fn, sites_fn, bc_fn){
  if(!file.size(mtx_fn) > 0){
    return(NULL)
  }

  rnames <- readLines(sites_fn)
  cnames <- readLines(bc_fn)
  dt <- data.table::fread(mtx_fn, sep = "\t", colClasses = "integer")

  if(ncol(dt) < 3) stop("malformed sparseMatrix")

  sp_idx <- 3:ncol(dt)
  lapply(sp_idx, function(x) {
    sparseMatrix(i = dt[[1]],
                 j = dt[[2]],
                 x = dt[[x]],
                 index1 = FALSE,
                 dims = c(length(rnames), length(cnames)),
                 dimnames = list(rnames, cnames))
  })
}

# Utilities -------------------------------------------------------

get_sc_pileup <- function(bamfn, id, sites, barcodes,
                          outfile_prefix, chrom,
                          umi_tag, cb_tag, libtype_code,
                          event_filters, fp, pe, verbose){
  if(length(chrom) > 0) {
    sites <- sites[seqnames(sites) %in% chrom, ]
  }

  if(length(sites) == 0) return(character())

  outfns <- c("counts.mtx", "sites.txt", "bcs.txt")
  plp_outfns <- file.path(outfile_prefix, paste0(id, "_", outfns))
  plp_outfns <- path.expand(plp_outfns)

  if(verbose) cli::cli_alert("working on group: {id}")
  lst <- gr_to_regions(sites)
  res <- .Call(".scpileup",
               bamfn,
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
               fp$min_variant_reads)

  if(res < 0) cli::cli_abort("pileup failed")
  plp_outfns
}

gr_to_regions <- function(gr){
  stopifnot(all(width(gr) == 1))

  nr <- length(gr);
  if(nr == 0)  cli::cli_abort("No entries in GRanges")

  gr_nms <- c("ref", "alt")
  if(!all(gr_nms %in% names(mcols(gr)))){
    cli::cli_abort("GRanges must have a ref and alt columns")
  }

  if(any(strand(gr) == "*")){
    cli::cli_alert_warning("missing strand not found in input, coercing strand to '+'")
    strand(gr) <- "+"
  }
  gr$idx <- seq(0, nr - 1) # is zero-based index


  list(as.character(seqnames(gr)),
       as.integer(start(gr)),
       as.integer(strand(gr)),
       as.character(mcols(gr)$ref),
       as.character(mcols(gr)$alt),
       as.integer(mcols(gr)$idx))

}

id_to_gr <- function(x, seq_info){
  xx <- str_split(x, ":", simplify = TRUE)
  if(ncol(xx) != 2) {
    cli::cli_abort("unable to decode sites to rownames")
  }
  seqnms <- xx[, 1]
  other_fields <- str_split(xx[, 2], "_", simplify = TRUE)

  if(ncol(other_fields) != 4) {
    stop("unable to decode sites to rownames")
  }

  gr <- GRanges(seqnms,
                IRanges(start = as.integer(other_fields[, 1]),
                        width = 1L),
                strand = c("+", "-")[as.integer(other_fields[, 2])],
                ref = other_fields[, 3],
                alt = other_fields[, 4],
                seqinfo = seq_info)
  names(gr) <- x
  gr
}

check_tag <- function(x) {
  if(length(x) != 1 && nchar(x) != 2){
    stop("supplied tag must by nchar of 2: ", x)
  }
}

chunk_vec <- function(x,n) {
  if(n == 1) {
    res <- list(`1` = x)
    return(res)
  }
  split(x, cut(seq_along(x), n, labels = FALSE))
}


