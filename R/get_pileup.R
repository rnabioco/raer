#' Generate pileup
#' @param bamfile path to bam file
#' @param fafile path to fasta file
#' @param bedfile path to bed file with sitesu or regions to query
#' @param region samtools region query string (i.e. chr1:100-1000)
#' @param chrom chrom to process, not to be used with region
#' @param start start position, not to be used with region
#' @param end end position, not to be used with region
#' @param library_type read orientation, one of fr-first-strand,
#' fr-second-strand, or unstranded.
#' @param outfile output tabix index file
#' @param return_data return data as a Granges table?
#'
#' @importFrom Rsamtools bgzip indexTabix TabixFile scanTabix
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @export
get_pileup <- function(bamfile,
                   fafile,
                   bedfile = NULL,
                   region = NULL,
                   chrom = NULL,
                   start = NULL,
                   end = NULL,
                   min_reads = 1L,
                   min_base_qual = 20L,
                   max_depth = 1e4,
                   library_type = "fr-first-strand",
                   outfile = NULL,
                   return_data = TRUE){

  bamfile <- path.expand(bamfile)
  fafile <- path.expand(fafile)

  if(!is.null(bedfile)){
    bedfile <- path.expand(bedfile)
    if(!file.exists(bedfile)){
      stop("bedfile not found: ", bedfile, call. = FALSE)
    }
  } else {
    bedfile = "."
  }

  if(!file.exists(bamfile)){
    stop("bamfile not found: ", bamfile, call. = FALSE)
  }

  if(!file.exists(fafile)){
    stop("fasta file not found: ", fafile, call. = FALSE)
  }

  if(is.null(outfile)){
    outfile <- tempfile()
  } else {
    outfile <- path.expand(outfile)
  }
  if(is.null(region)) {
    if(!is.null(chrom)){
      region = chrom
      if(!is.null(start) & !is.null(end)){
        region = paste0(region, ":", start, end)
      }

    } else {
      # not sure how to catch NULL in rcpp so using "." to indicate no region instead
      region = "."
    }
  }

  lib_types <- c("fr-first-strand", "fr-second-strand", "unstranded")
  if(!library_type %in% lib_types){
    stop("library_type must be one of :", paste(lib_types, collapse = " "))
  }

  res <- run_pileup(bamfile,
                    fafile,
                    region,
                    outfile,
                    bedfile,
                    min_reads,
                    max_depth,
                    min_base_qual,
                    library_type)

  if(file.info(outfile)$size == 0){
    return()
  }

  tbxfile <- Rsamtools::bgzip(outfile, overwrite = TRUE)
  idx <- Rsamtools::indexTabix(tbxfile, seq = 1, start = 2, end = 2, zeroBased = FALSE)
  tbx <- Rsamtools::TabixFile(tbxfile, idx)

  if(region != "."){
    ivl_vals <- get_region(region)
    params = GenomicRanges::GRanges(ivl_vals$chrom,
                                    IRanges::IRanges(start = ivl_vals$start + 1,
                                                     end = min(ivl_vals$end, MAX_INT)))
    tbx_vals <- Rsamtools::scanTabix(tbx, param = params)[[1]]
  } else {
    tbx_vals <- Rsamtools::scanTabix(tbx)[[1]]
  }
  from <- data.table::fread(paste0(paste(tbx_vals,collapse="\n"),"\n" ),
                stringsAsFactors = FALSE,
                data.table = FALSE,
                sep="\t")

  colnames(from)[4:ncol(from)] = c("Ref", "nRef", "nVar",
                                   "nA", "nT", "nC",
                                   "nG", "nN")
  if(return_data){
    GenomicRanges::GRanges(seqnames=from$V1,
                           ranges=IRanges::IRanges(start=from$V2,
                                          end=from$V2 + 1),
                           strand=from$V3,
                           from[,4:ncol(from)])
  }
}

####################### Utils

MAX_INT <- 536870912

