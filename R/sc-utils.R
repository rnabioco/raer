
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
                         pos_sort_output = TRUE){
  stopifnot(file.exists(inbam))
  inbam <- path.expand(inbam)
  idx_file <- paste0(inbam, ".bri")
  if(!file.exists(idx_file)){
    stop("bam file must be sorted by tag, and indexed with build_index\n",
         "samtools sort -t CB your.bam\n",
         "then index in R:\n",
         "raer::build_tag_index('your_sorted.bam', tag = 'CB')")
  }
  stopifnot(is.character(barcodes))

  tmp_bam <- tempfile(fileext = ".bam")
  if(is.null(outbam)){
    outbam <- tempfile(fileext = ".bam")

  }

  fetch_cb_reads(inbam, tmp_bam, barcodes);

  if(pos_sort_output){
    # Rsamtools will add .bam to end of the output bam file
    outbam <- gsub(pattern = ".bam$", "", outbam)
    outbam <- Rsamtools::sortBam(tmp_bam, outbam)
    Rsamtools::indexBam(outbam)
  }

  unlink(tmp_bam)
  return(outbam)
}
