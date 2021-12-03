#' Generate pileup
#' @param bamfile path to bam file
#' @param fafile path to fasta file
#' @param bedfile path to bed file with sites to query
#' @param region samtools region query string (i.e. chr1:100-1000)
#' @param chrom chrom to process, not to be used with region
#' @param start start position, not to be used with region
#' @param end end position , not to be used with region
#' @param outfile output tabix index file
#' @param return_data return data as a Granges table?
#'
#' @importFrom Rsamtools bgzip indexTabix TabixFile scanTabix
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @export
get_pileup <- function(bamfile = NULL,
                   fafile = NULL,
                   bedfile = NULL,
                   region = NULL,
                   chrom = NULL,
                   start = NULL,
                   end = NULL,
                   outfile = "tmp.txt",
                   return_data = TRUE){
  if(is.null(bamfile)) {
    bamfile <- system.file("extdata",
                           "SRR1258218_chr21.bam",
                           package = "rino")
  }

  if(is.null(fafile)) {
    fafile <- system.file("extdata",
                          "chr21.fa",
                          package = "rino")
  }

  if(is.null(bedfile)) {
    bedfile <- system.file("extdata",
                          "chr21_regions.bed",
                          package = "rino")
  }

  bamfile <- path.expand(bamfile)
  fafile <- path.expand(fafile)
  outfile <- path.expand(outfile)
  bedfile <- path.expand(bedfile)

  if(is.null(region)) {
    if(any(!c(is.null(chrom), is.null(start), is.null(end)))){
      region = paste0(chrom, ":", start, end)
    } else {
      # not sure how to pass NULL to rcpp so using "." to indicate no region instead
      region = "."
    }
  }

  res <- run_pileup(bamfile, fafile, region, outfile, bedfile)

  tbxfile <- Rsamtools::bgzip(outfile, overwrite = TRUE)
  Rsamtools::indexTabix(tbxfile, seq = 1, start = 2, end = 2, zeroBased = FALSE)
  tbx <- Rsamtools::TabixFile(tbxfile)

  if(region != "."){
    ivl_vals <- get_region(region)
    params = GenomicRanges::GRanges(ivl_vals$chrom,
                                    IRanges::IRanges(start=ivl_vals$start + 1,
                                                     end=min(ivl_vals$end, 536870912)))
    tbx_vals <- Rsamtools::scanTabix(tbx, param = params)[[1]]
  } else {
    tbx_vals <- Rsamtools::scanTabix(tbx)[[1]]
  }
  from <- data.table::fread(paste0(paste(tbx_vals,collapse="\n"),"\n" ),
                stringsAsFactors = FALSE,
                data.table = FALSE,
                sep="\t")

  if(return_data){
    GenomicRanges::GRanges(seqnames=from$V1,
                           ranges=IRanges::IRanges(start=from$V2,
                                          end=from$V2 + 1),
                           strand=from$V3,
                           from[,4:ncol(from)])
  }
}
