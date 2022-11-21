#' A reference class for generating and storing index of BedFile
#'
#' @description Class to store information about a bedfile used
#' for indicating regions for pileup
#'
#' @field .extptr an externalptr to c-level samtools bedindex
#' @field path filepath to bed file
#' @field open logical indicating if the index is open
#' @rdname index_bed
.BedFile <- setRefClass("BedFile",
  fields = list(
    .extptr = "externalptr",
    path = "character",
    open = "logical"
  )
)
