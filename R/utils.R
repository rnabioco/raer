
is_null_extptr <- function(pointer) {
  stopifnot(is(pointer, "externalptr"))
  .Call(".isnull", pointer)
}


#' Read in tabix indexed file as a data.frame
#'
#' @param filename path to tabix file
#' @param region samtools region query string (i.e. chr1:100-1000)
#' @param numeric_cols columns to convert to numeric as a integer vector (one-based index)
#' @param col_names column names in the output data.frame

#' @examples
#' bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#'
#' tbx_fn <- get_pileup(bamfn, fafn, return_data = FALSE)
#' head(read_tabix(tbx_fn))
#' read_tabix(tbx_fn, region = "SPCS3:498-500")
#' @rdname read_tabix
#' @export
read_tabix <- function(filename,
                       region = ".",
                       numeric_cols = c(2, 6:12),
                       col_names = PILEUP_COLS){
  filename <- path.expand(filename)
  # returned as a list to avoid stringsAsFactors
  df <- cread_tabix(filename, region)
  numeric_cols <- intersect(c("start", "end", "pos"),  colnames(df))
  if(!is.null(numeric_cols)){
    for(i in numeric_cols){
      df[[i]] <- as.numeric(df[[i]])
    }
  }
  colnames(df) <- col_names

  df
}

PILEUP_COLS <-  c("chrom",
                  "pos",
                  "strand",
                  "Ref",
                  "Var",
                  "nRef",
                  "nVar",
                  "nA",
                  "nT",
                  "nC",
                  "nG",
                  "nN")

#' List chromosomes in a tabix index
#' @param filename path to indexed tabix file
#'
#' @examples
#' bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#'
#' tbx_fn <- get_pileup(bamfn, fafn, return_data = FALSE)
#' get_tabix_chroms(tbx_fn)
#' @rdname read_tabix
#' @export
get_tabix_chroms <- function(filename){
  list_tabix_chroms(path.expand(filename))
}

#' Provide working directory for raer example files.
#'
#' @param path path to file
#'
#' @examples
#' raer_example('human.fasta')
#'
#' @export
raer_example <- function(path) {
  system.file("extdata", path, package = "raer", mustWork = TRUE)
}


#' Copy of VariantTools::extractCoverageForPositions
#' Authors: Michael Lawrence, Jeremiah Degenhardt, Robert Gentleman
#' @param cov coverage produced by [GenomicAlignments::coverage()]
#' @param pos GRanges containing editing sites
#' @importFrom GenomeInfoDb `seqlevels<-` seqlevels
getCoverageAtPositions <- function(cov, pos){
  if (length(setdiff(seqlevels(pos), names(cov))) > 0L)
    stop("Some seqlevels are missing from coverage")
  if (any(width(pos) > 1L))
    stop("Some ranges are of width > 1")
  seqlevels(pos) <- names(cov)
  ord <- order(seqnames(pos))
  ans <- integer(length(pos))
  ans[ord] <- unlist(mapply(function(v, p) {
    runValue(v)[findRun(p, v)]
  }, cov, split(start(pos), seqnames(pos)), SIMPLIFY=FALSE), use.names=FALSE)
  ans
}
