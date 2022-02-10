
#' Read in bam file as a data.frame
#' @param filename path to bam file
#' @param region samtools region query string (i.e. chr1:100-1000)
#' @param tags bam tags to return, supplied with type suffix, i.e. "XS:A" for character,
#' "AS:i" for integer, and "CB:Z" for string. Bam entries with missing tags will return empty
#' strings.
#' @examples
#' bam_file <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")
#' res <- bam_to_df(bam_file)
#' head(res)
#'
#' res <- bam_to_df(bam_file,
#'                  region = "DHFR:88-100",
#'                  tags = c("NH:i", "AS:i"))
#' res
#' @export
bam_to_df <- function(filename,
                      region = NULL,
                      tags = NULL){

  filename <- path.expand(filename)
  if(!file.exists(filename)){
    stop("file not found: ", filename)
  }
  if(!is.null(tags)){
    #parse tag type
    if(any(!grepl(":[ZiA]$", tags))){
      stop("missing type value in supplied tags, need :Z, :i, or :A suffix on tags")
    }
    sp_tags <- strsplit(tags, ":")
    tag_types <- unlist(lapply(sp_tags, function(x) x[2]))
    tag_ids <- unlist(lapply(sp_tags, function(x) x[1]))

  } else{
    tag_types = ""
    tag_ids = ""
  }

  if(is.null(region)){
    region = "."
  }
  res <- read_bam(filename, tag_ids, tag_types, region)

  res
}
