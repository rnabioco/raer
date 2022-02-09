
#' Subset a bam file to contains only certain cell barcodes
#' @export
get_cell_bam <- function(inbam,
                         barcodes = c("AAACCTGGTGACAAAT-1"),
                         outbam = NULL){
  stopifnot(file.exists(inbam))
  idx_file <- paste0(inbam, ".bri")
  if(!file.exists(idx_file)){
    stop("bam file must be sorted by cellbarcode, and indexed with build_index")
  }
  stopifnot(is.character(barcodes))

  tmp_bam <- tempfile(fileext = ".bam")
  if(is.null(outbam)){
    outbam <- tempfile(fileext = ".bam")

  }
  message(tmp_bam)
  message(outbam)
  raer:::fetch_cb_reads(inbam,  tmp_bam,   barcodes);
  Rsamtools::sortBam(tmp_bam, outbam)
  Rsamtools::indexBam(outbam)
  unlink(tmp_bam)
  return(outbam)
}
