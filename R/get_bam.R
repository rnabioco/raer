
#' Read in bam file as a data.frame
#' @param filename path to bam file
#' @param region samtools region query string (i.e. chr1:100-1000)
#' @param tags bam tags to return, supplied with type suffix, i.e. "XS:A" for character,
#' "AS:i" for integer, and "CB:Z" for string. Bam entries with missing tags will return empty
#' strings.
#' @examples
#' bam_file <- system.file("extdata", "small_sorted.bam", package = "kentr")
#' res <- bam_to_df(bam_file)
#' head(res)
#'
#' res <- bam_to_df(bam_file,
#'                  region = "chr1:3e6-3.1e6",
#'                  tags = c("XS:A", "AS:i", "CB:Z"))
#' res
#' @export
bam_to_df <- function(filename = NULL,
                      region = ".",
                      tags = NULL){
  if(is.null(filename)) filename <- system.file("extdata",
                                                "small_sorted.bam",
                                                package = "rino")

  filename <- path.expand(filename)


  if(!is.null(tags)){
    #parse tag type
    if(any(!str_detect(tags, ":[ZiA]$"))){
      stop("missing type value in supplied tags, need :Z, :i, or :A suffix on tags")
    }
    sp_tags <- stringr::str_split(tags, ":", simplify = T)
    tag_types <- sp_tags[, 2]
    tag_ids <- sp_tags[, 1]

  } else{
    tag_types = ""
    tag_ids = ""
  }
  res <- read_bam(filename, region, tag_ids, tag_types)

  res
}
