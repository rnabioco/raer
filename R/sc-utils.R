#' Build index for tag sorted bam file
#'
#' @param bamfile tag sorted bamfile
#' @param tag name of tag in bamfile used for sorting. The tag must be of type
#'   "Z".
#' @param n_records_to_check The number of bam records to query to validate that
#'   the tag is present and of the correct type. Set to 0 to disable checks.
#' @param overwrite if TRUE, regenerate index if it already exists
#' @return name of index generated, which is the bam file + ".bri"
#'
#' @examples
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
#' build_tag_index(bam_fn)
#'
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_ubsort.bam", package = "raer")
#' build_tag_index(bam_fn, tag = "UB")
#'
#' @importFrom Rsamtools scanBamHeader scanBam isOpen BamFile yieldSize<-
#'   ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#' @export
build_tag_index <- function(bamfile, tag = "CB", n_records_to_check = 1e6,
                            overwrite = TRUE) {
  bamfile <- path.expand(bamfile)
  stopifnot(file.exists(bamfile))
  stopifnot(is.character(tag) && length(tag) == 1 && nchar(tag) == 2)

  outfn <- paste0(bamfile, ".bri")
  if (!overwrite && file.exists(outfn)) {
    return(outfn)
  }

  bfo <- Rsamtools::BamFile(bamfile)
  so <- Rsamtools::scanBamHeader(bfo)$text$`@HD`[2]
  if (so == "SO:coordinate") {
    stop("bam file must be sorted by tag, not coordinate sorted")
  }

  if (n_records_to_check < 1) {
    message("disabling checks for tag in bam\n",
      "set n_records_to_check > 0 to enable checks")
  } else {
    Rsamtools::yieldSize(bfo) <- n_records_to_check
    mapped <- scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE)
    alns <- Rsamtools::scanBam(bfo,
      param = Rsamtools::ScanBamParam(tag = tag,
        flag = mapped))
    tag_vals <- alns[[1]]$tag[[tag]]
    if (length(tag_vals) == 0) {
      stop("unable to find tag ", tag,
        " in first ", n_records_to_check,
        " records of bam file")
    }
    if (!is.character(tag_vals)) {
      stop("invalid type ", typeof(tag_vals), " \n",
        "Indexing only supports type character\n",
        "e.g: tags with Z, e.g. CB:Z:ABCBACB")
    }
  }
  if (Rsamtools::isOpen(bfo)) close(bfo)

  c_build_index(bamfile, outfn, tag)
  return(outfn)
}

#' Show tags stored in tag index
#'
#' @param bamfile tag sorted bamfile, indexed with [build_tag_index()]
#' @return Character vector of tags
#'
#' @examples
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
#' build_tag_index(bam_fn)
#' show_tag_index(bam_fn)
#'
#' @export
show_tag_index <- function(bamfile) {
  bamfile <- path.expand(bamfile)
  idxfn <- paste0(bamfile, ".bri")
  stopifnot(file.exists(bamfile))
  stopifnot(file.exists(idxfn))
  c_show_index(bamfile, idxfn)
}


#' Subset a bam file to contain only certain tags (e.g. cell barcodes)
#'
#' @param bamfile input tag indexed bam file
#' @param barcodes character vector of tag values to extract
#' @param outbam optional output bam file name, if not supplied a temporary file
#'   will be used
#' @param pos_sort_output if TRUE, sort output bamfile by position and generate
#'   a samtools style index
#' @param ... Additional arguments passed to [Rsamtools::sortBam()]
#'
#' @returns Returns name of output bam file. the output bam file will be
#'   positionally sorted and positionally indexed using Rsamtools.
#'
#' @examples
#' library(GenomicAlignments)
#'
#' bam_fn <- raer_example("5k_neuron_mouse_xf25_1pct_cbsort.bam")
#' build_tag_index(bam_fn)
#'
#' cbs <- c("AGGATAATCTCAGAAC-1", "TTCGATTTCCCGAGGT-1")
#' bam_out <- get_tag_bam(bam_fn, barcodes = cbs)
#' readGAlignments(bam_out, param = ScanBamParam(tag = "CB"))
#'
#' @export
get_tag_bam <- function(bamfile,
                        barcodes,
                        outbam = NULL,
                        pos_sort_output = TRUE,
                        ...) {
  stopifnot(file.exists(bamfile))
  bamfile <- path.expand(bamfile)
  idx_file <- paste0(bamfile, ".bri")
  if (!file.exists(idx_file)) {
    stop("bam file must be sorted by tag, and indexed with build_index\n",
      "samtools sort -t CB your.bam\n",
      "then index in R:\n",
      "raer::build_tag_index('your_sorted.bam', tag = 'CB')")
  }
  stopifnot(is.character(barcodes) && length(barcodes) > 0)

  tmpbam <- tempfile(fileext = ".bam")
  if (is.null(outbam)) {
    outbam <- tempfile(fileext = ".bam")
  }

  fetch_cb_reads(bamfile, tmpbam, barcodes)

  if (pos_sort_output) {
    # tag sorted bams are secondarily sorted by position
    # so no need to re-sort if grabbing a single tag
    if (length(barcodes) > 1) {
      # Rsamtools will add .bam to end of the output bam file
      outbam <- gsub(pattern = ".bam$", "", outbam)
      outbam <- Rsamtools::sortBam(tmpbam, outbam, ...)
    } else {
      file.rename(tmpbam, outbam)
    }
    Rsamtools::indexBam(outbam)
  } else {
    file.rename(tmpbam, outbam)
  }

  unlink(tmpbam)
  return(outbam)
}


# validate tags are in index
check_missing_barcodes <- function(cbs, bamfile) {
  tags <- show_tag_index(bamfile)
  cbs <- unlist(cbs, use.names = FALSE)
  sum(!cbs %in% tags$tag)
}

#' Identify sites with read coverage
#'
#' @description This function will compute read coverage at sites
#' listed in either a GRanges or BED file. Sites passing with at least
#' `min_counts` read depth will be returned. This function will
#' compute coverage across all alignments and positions,
#' ignoring any cell-barcodes. Generally it is faster to compute coverage
#' across all positions, then filter to those in the GRanges object, rather
#' than only querying the supplied sites.
#'
#' @param bamfile Bam file name
#' @param gr GRanges, or file name of BED file, containing site to
#' examine
#' @param min_counts minimum reads counts to require
#' @param param A ScanBamParam object specifying how reads should be filtered.
#' If not supplied the default behavior will ignore exclude alignments marked
#' as secondary, duplicate, supplementary, or QC-fail.
#' @param verbose If TRUE, print messages
#' @param ... Additional arguments to supply to [GenomicAlignments::coverage()]
#'
#' @importFrom GenomicAlignments coverage
#' @export
filter_by_coverage <- function(bamfile, gr, min_counts,
                               param = NULL, verbose = FALSE,
                               ...) {
  if(is.character(gr)){
    if(!file.exists(gr)) stop(gr, " file not found")
    gr <- rtracklayer::import(gr)
  }
  gr <- GenomicRanges::sort(gr)
  n_sites <- length(gr)

  if(is.null(param)){
    covflags <- scanBamFlag(isSecondaryAlignment = FALSE,
                isDuplicate = FALSE,
                isSupplementaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE)
    param <- ScanBamParam(flag = covflags)
  }

  cvg <- GenomicAlignments::coverage(bamfile, param = param, ...)
  shared_seqs <- intersect(names(cvg), GenomeInfoDb::seqlevels(gr))
  if(length(shared_seqs) == 0){
    stop("no shared seqnames found, check input bamfile and bedfile chromosome names")
  }
  cvg <- cvg[shared_seqs]
  gr <- gr[seqnames(gr) %in% shared_seqs, ]
  seqlevels(cvg) <- shared_seqs
  seqlevels(gr) <- shared_seqs
  cvg <- getCoverageAtPositions(cvg, gr)
  mcols(gr)$score <- cvg
  gr <- gr[cvg >= min_counts]
  if (verbose) {
    message("Input bed contained ", n_sites, " sites\n",
            n_sites - length(gr), " sites do not have at least ",
            min_counts, " reads \n",
            "and will be ignored.")
  }
  gr
}


get_cell_pileup <- function(bamfn, fafn, cellbarcodes, per_cell, ...) {
  if(per_cell){
    bam <- lapply(cellbarcodes,
                          function(x){
                            get_tag_bam(bamfn,
                                        barcodes = x,
                                        outbam = NULL)})
    bam <- unlist(bam)
  } else {
    bam <- get_tag_bam(bamfn,
                       barcodes = cellbarcodes,
                       outbam = NULL)
  }
  on.exit(unlink(c(bam, paste0(bam, ".bai"))))
  out <- get_pileup(bam, fafile = fafn, ...)
  out
}


#' Calculate editing frequencies per cell or per cluster
#'
#' @param bamfile BAM file name
#' @param fafile FASTA file name
#' @param bedfile BED file containing editing sites
#' @param cell_barcodes A list of character vectors containing cell barcodes
#'   to query or a single character vector of single cells to process. If a list
#'   is supplied it is assumed that alignments from cell barcodes in each list element
#'   should be pooled. If a single character vector, it is assumed that each cell barcode
#'   should be processed independently. See examples for specification.
#' @param min_reads Minimum read counts required to consider a site for editing.
#'   This is calculated across all reads in the bamfile prior to running `[get_pileup()]`
#'   per cell to remove low-frequency events. If set to 0 (default) this step will
#'   not be run. See `filter_by_coverage()` for more details.
#' @param assay_cols assays to store in returned se. Set to "A" and "G". Note that
#'   storing multiple assays can require large amounts of memory.
#' @param tag_index_args arguments pass to [`build_tag_index()`]
#' @param sparse if TRUE, store matrices in sparseMatrix format in SummarizedExperiment.
#' @param BPPARAM BiocParallel instance. Parallel computation occurs across
#'   each entry in the cell_barcodes list, or across batches of single cells specified
#'   by batch_size.
#' @param batch_size When processing single cells, the batch_size controls
#'   how many individual cell bams to process in each invocation of `get_pileup()`.
#'   Batching the cells reduces run time by avoiding  loading sequences from the
#'   fasta file for each cell. Setting values above 50 is unlikely to further improve
#'   runtime.
#' @param verbose Display messages
#' @param ... additional arguments passed to `[get_pileup()]`.
#'
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' # get vector of cell barcodes in bam file (for use in this example)
#' # usually these would come from the single cell analysis
#'
#' bamfn <- raer_example("5k_neuron_mouse_xf25_1pct_cbsort.bam")
#' idxfn <- build_tag_index(bamfn)
#' cbs <- show_tag_index(bamfn)$tag
#' cbs[1:5]
#'
#' # process each cell individually
#' # will be slow with many sites and cells
#' # bam file will be indexed by build_tag_index() if not already done.
#' fp <- FilterParam(library_type = "fr-second-strand")
#' se <- sc_editing(bamfile = bamfn,
#'   fafile = raer_example("mouse_tiny.fasta"),
#'   bedfile = raer_example("5k_neuron_sites.bed.gz"),
#'   cell_barcodes = cbs[1:15],
#'   verbose = FALSE,
#'   filterParam = fp)
#'
#' # pool cell barcodes across clusters
#' # pass a named list, with each list entry corresponding to a vector
#' # of cell barcodes from a cluster
#'
#' # simulate 5 clusters
#' cb_lst <- split(cbs, cut(seq_along(cbs), breaks = 5))
#' names(cb_lst) <- paste0("cluster", 1:5)
#'
#' se <- sc_editing(bamfile = bamfn,
#'   fafile = raer_example("mouse_tiny.fasta"),
#'   bedfile = raer_example("5k_neuron_sites.bed.gz"),
#'   cell_barcodes = cb_lst,
#'   verbose = FALSE,
#'   filterParam = fp)
#' assays(se)$nA
#' assays(se)$nG
#'
#' @importFrom rtracklayer import
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom methods formalArgs

#' @export
sc_editing <- function(bamfile,
                       fafile,
                       bedfile,
                       cell_barcodes,
                       min_reads = 0L,
                       assay_cols = c("nA", "nG"),
                       tag_index_args = list(tag = "CB"),
                       sparse = TRUE,
                       BPPARAM = SerialParam(),
                       batch_size = 50,
                       verbose = TRUE,
                       ...) {

  if (!all(assay_cols %in% PILEUP_COLS)) {
    allowed_vals <- paste(PILEUP_COLS[5:length(PILEUP_COLS)],
      collapse = ", ")
    stop("assay_cols input not correct\n  ",
      "must match ", allowed_vals)
  }

  # fail early if incorrect args passed through ...
  plp_args <- names(list(...))
  invalid_args <- plp_args[!plp_args %in% formalArgs(get_pileup)]
  if (length(invalid_args) > 0) {
    stop(invalid_args, " is not a valid argument for get_pileup()")
  }
  idx_fn <- paste0(bamfile, ".bri")
  if (!file.exists(idx_fn)) {
    if (verbose) message("building cellbarcode index for bam file")
    do.call(build_tag_index, c(bamfile = bamfile, tag_index_args))
  }

  if (!is.list(cell_barcodes)) {
    per_cell <- TRUE
    stopifnot(batch_size > 0)
    ids <- cell_barcodes
    cell_barcodes <- split_vec(cell_barcodes, batch_size)
    #names(cell_barcodes) <- ids
  } else {
    per_cell <- FALSE
  }

  n_invalid_bcs <- check_missing_barcodes(cell_barcodes, bamfile)
  if (n_invalid_bcs > 0) {
    warning(n_invalid_bcs, " cell_barcodes are missing from the tag index.")
  }

  if (min_reads > 0) {
    if (verbose) message("Examining coverage at supplied sites.")
    bed <- filter_by_coverage(bamfile, bedfile, min_reads, verbose = verbose)
    if (length(bed) == 0) stop("no sites remaining to process")
    tmp_bed <- tempfile(fileext = ".bed")
    on.exit(unlink(tmp_bed), add = TRUE)
    export(bed, tmp_bed)
    bedfile <- tmp_bed
  }

  # build c-level hash of regions to keep once, rather at every iteration
  # probably won't work on windows.
  idx <- indexBed(bedfile)
  on.exit(close(idx), add = TRUE)

  if (verbose) message("beginning pileup")
  res <- bplapply(seq_along(cell_barcodes), function(i) {
    if (verbose) {
      if(per_cell){
        message("working on: batch ", i,
                " cells [", 1 + ((i-1) * batch_size),
                       "-",
                i * batch_size,
                      "]")
      } else {
        message("working on: group ", i)
      }
    }
    plps <- get_cell_pileup(bamfile, fafile, cell_barcodes[[i]],
                    per_cell = per_cell,
      bedfile = NULL, bedidx = idx, return_data = TRUE,
      verbose = verbose,
      BPPARAM = BiocParallel::SerialParam(),
      ...)
    if(per_cell){
      names(plps) <- cell_barcodes[[i]]
    }
    plps
  }, BPPARAM = BPPARAM)
  bpstop(BPPARAM)

  if(per_cell) {
    res <- unlist_w_names(res)
    # figure out numeric assays
    num_cols <- assay_cols[which(sapply(mcols(res[[1]])[assay_cols], is.numeric))]
    # remove zero depth (helps with making the sparseMatrices)
    res <- lapply(res, function(x) {
      to_keep <- rowSums(as.matrix(mcols(x)[num_cols])) > 0
      x[to_keep]
    })
   # res <- lapply(res, function(x), )
  } else {
    names(res) <- names(cell_barcodes)
  }

  if (verbose) message("collecting pileups into summarizedExperiment")
  se <- create_se(res,
    assay_cols = assay_cols,
    sparse = sparse,
    fill_na = 0L,
    verbose = verbose)
  se
}
