
#' Build index for tag sorted bam file
#' @param bamfn tag sorted bamfile
#' @param tag name of tag in bamfile used for sorting. The tag must be of type "Z".
#' @param n_records_to_check The number of bam records to query to validate that the
#' tag is present and of the correct type. Set to 0 to disable checks.
#' @return name of index generated, which is the bam file + ".bri"
#'
#' @examples
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
#' build_tag_index(bam_fn)
#'
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_ubsort.bam", package = "raer")
#' build_tag_index(bam_fn, tag = "UB")
#'
#' @importFrom Rsamtools scanBamHeader scanBam isOpen BamFile yieldSize<- ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#' @export
build_tag_index <- function(bamfn, tag = "CB", n_records_to_check = 1e3){
  stopifnot(file.exists(bamfn))
  stopifnot(is.character(tag) && length(tag) == 1 && nchar(tag) == 2)

  bamfn <- path.expand(bamfn)
  outfn <- paste0(bamfn, ".bri")

  bfo <- Rsamtools::BamFile(bamfn)
  so <- Rsamtools::scanBamHeader(bfo)$text$`@HD`[2]
  if(so == "SO:coordinate"){
    stop("bam file must be sorted by tag, not coordinate sorted")
  }

  if(n_records_to_check < 1){
    message("disabling checks for tag in bam\n",
            "set n_records_to_check > 0 to enable checks")
  } else {
    Rsamtools::yieldSize(bfo) <- n_records_to_check
    alns <- Rsamtools::scanBam(bfo, param = Rsamtools::ScanBamParam(tag = tag))
    tag_vals <- alns[[1]]$tag[[tag]]
    if(length(tag_vals) == 0){
      stop("unable to find tag ", tag,
           " in first ", n_records_to_check,
           " records of bam file")
    }
    if(!is.character(tag_vals)){
      stop("invalid type ", typeof(tag_vals), " \n",
           "Indexing only supports type character\n",
           "e.g: tags with Z, e.g. CB:Z:ABCBACB")
    }
  }
  if(Rsamtools::isOpen(bfo)) close(bfo)

  c_build_index(bamfn, outfn, tag)
  return(outfn)
}


#' Subset a bam file to contain only certain tags (e.g. cell barcords)
#'
#' @param inbam input tag indexed bam file
#' @param barcodes character vector of tag values to extract
#' @param outbam optional output bam file name, if not supplied a temporary
#' file will be used
#' @param pos_sort_output if TRUE, sort output bamfile by position and generate
#' a samtools style index
#' @param ... Additional arguments passed to [Rsamtools::sortBam()]
#'
#' @returns Returns name of output bam file. the output bam file will be
#' positionally sorted and positionally indexed using Rsamtools.
#'
#' @examples
#' library(GenomicAlignments)
#' bam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
#' build_tag_index(bam_fn)
#' cbs <- c("AGGATAATCTCAGAAC-1", "TTCGATTTCCCGAGGT-1")
#' bam_out <- get_cell_bam(bam_fn, barcodes = cbs)
#' readGAlignments(bam_out, param = ScanBamParam(tag = "CB"))
#' @export
get_cell_bam <- function(inbam,
                         barcodes,
                         outbam = NULL,
                         pos_sort_output = TRUE,
                         ...){
  stopifnot(file.exists(inbam))
  inbam <- path.expand(inbam)
  idx_file <- paste0(inbam, ".bri")
  if(!file.exists(idx_file)){
    stop("bam file must be sorted by tag, and indexed with build_index\n",
         "samtools sort -t CB your.bam\n",
         "then index in R:\n",
         "raer::build_tag_index('your_sorted.bam', tag = 'CB')")
  }
  stopifnot(is.character(barcodes) && length(barcodes) > 0)

  tmp_bam <- tempfile(fileext = ".bam")
  if(is.null(outbam)){
    outbam <- tempfile(fileext = ".bam")

  }

  fetch_cb_reads(inbam, tmp_bam, barcodes);

  if(pos_sort_output){
    # Rsamtools will add .bam to end of the output bam file
    outbam <- gsub(pattern = ".bam$", "", outbam)
    outbam <- Rsamtools::sortBam(tmp_bam, outbam, ...)
    Rsamtools::indexBam(outbam)
  }

  unlink(tmp_bam)
  return(outbam)
}

#' @importFrom GenomicAlignments coverage
#' @importFrom VariantTools extractCoverageForPositions
filter_by_coverage <- function(bamfile, gr, min_counts, ...){
  cov <- GenomicAlignments::coverage(bamfile, ...)
  cov <- cov[names(cov) %in% seqlevels(gr)]
  cov <- VariantTools::extractCoverageForPositions(cov, gr)
  gr[cov >= min_counts]
}


get_cell_pileup <- function(bamfn, fafn, cellbarcodes, ...){
  cluster_bam <- get_cell_bam(bamfn,
                              barcodes = cellbarcodes,
                              NULL,
                              maxMemory = 1024)
  on.exit(unlink(c(cluster_bam, paste0(cluster_bam, ".bri"))))

  out <- get_pileup(cluster_bam,
                    fafile = fafn,
                    bedfile = NULL,
                    ...)

  unlink(cluster_bam)
  out
}


#' Calculate editing frequencies
#'
#' @param bamfn BAM file name
#' @param fafn FASTA file name
#' @param bedfn BED file containing editing sites
#' @param cell_barcodes List of character vectors containing cell barcodes
#' to query. See examples for specification.
#' @param cell_bc_tag tag in bam file containing cell barcodes
#' @param min_counts Minimum read counts required to consider a site for editing.
#' This is calculated across all reads in the bamfile prior to running `get_pileup()`
#' per cell to remove low-frequency events. Read alignments are required to be
#' non-duplicate, primary, QC passing, and non-supplemental.
#' @param assay_cols assays to store in returned se. Set to "A" and "G". Note that
#' storing multiple assays can require large amounts of memory.
#' @param ... additional arguments passed to `[get_pileup()]`.
#' @param BPPARAM BiocParallel instance. Parallel computation occurs across
#' each entry in the cell_barcodes list
#' @param verbose Display messages
#'
#'
#' @importFrom rtracklayer import
#' @importFrom Rsamtools ScanBamParam scanBamFlag

#' @export
sc_editing <- function(bamfn,
                       fafn,
                       bedfn,
                       cell_barcodes,
                       cell_bc_tag = "CB",
                       min_counts = 25L,
                       assay_cols = c("nA", "nG"),
                       ...,
                       BPPARAM = SerialParam(),
                       verbose = TRUE){
  if(min_counts > 0){
    if(verbose) message("calculating coverage using all reads")
    bed <- rtracklayer::import(bedfn)
    n_sites <- length(bed)
    bed <- filter_by_coverage(bamfn, bed, min_counts,
                              param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                                                      isDuplicate = FALSE,
                                                                      isSupplementaryAlignment = FALSE,
                                                                      isNotPassingQualityControls = FALSE)))

    message("Input bed contained ", n_sites, " sites\n",
            "After filtering by min_count, ", length(bed), " sites remain")
    tmp_bed <- tempfile(fileext = ".bed")
    on.exit(unlink(tmp_bed))
    export(bed, tmp_bed)
    bedfn <- tmp_bed
  }

  if(!file.exists(paste0(bamfn, ".bri"))){
    if(verbose) message("building cellbarcode index for bam file")
    build_tag_index(bamfn)
  }

  # build c-level hash of regions to keep once, rather at every iteration
  idx <- indexBed(bedfn)

  if(verbose) message("beginning pileup")
  res <- bplapply(seq_along(cell_barcodes), function(i){
    if(verbose){
      message("working on: ", names(cell_barcodes)[i])
    }
    get_cell_pileup(bamfn, fafn, cell_barcodes[[i]],
                    bedidx = idx, return_data = TRUE, ...)
  }, BPPARAM = BPPARAM)
  bpstop(BPPARAM)

  if(verbose) message("collecting pileups into summarizedExperiment")
  se <- create_se(res,
                  assay_cols = assay_cols,
                  sparse = TRUE,
                  fill_na = 0L)
  se
}
