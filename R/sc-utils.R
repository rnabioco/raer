
#' Build index for tag sorted bam file
#' @param bamfn CB tag sorted bamfile
#' @return name of index generated, which is the bam file + ".bri"
#' @export
build_tag_index <- function(bamfn){
  stopifnot(file.exists(bamfn))
  bamfn <- path.expand(bamfn)
  outfn <- paste0(bamfn, ".bri")
  c_build_index(bamfn, outfn)
  return(outfn)
}


#' Subset a bam file to contain only certain cell barcodes
#' @param inbam input tag indexed bam file
#' @param barcodes character vector of tag values to extract
#' @param outbam optional output bam file name
#'
#' @return  Returns name of output bam file
#' @export
get_cell_bam <- function(inbam,
                         barcodes,
                         outbam = NULL){
  stopifnot(file.exists(inbam))
  inbam <- path.expand(inbam)
  idx_file <- paste0(inbam, ".bri")
  if(!file.exists(idx_file)){
    stop("bam file must be sorted by CB tag, and indexed with build_index")
  }
  stopifnot(is.character(barcodes))

  tmp_bam <- tempfile(fileext = ".bam")
  if(is.null(outbam)){
    outbam <- tempfile(fileext = ".bam")

  }

  fetch_cb_reads(inbam, tmp_bam, barcodes);

  # Rsamtools will add .bam to end of the output bam file
  outbam <- gsub(pattern = ".bam$", "", outbam)
  outbam <- Rsamtools::sortBam(tmp_bam, outbam)
  Rsamtools::indexBam(outbam)
  unlink(tmp_bam)
  return(outbam)
}
