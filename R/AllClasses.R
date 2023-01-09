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

#' A reference class for generating and storing region index
#'
#' @description Class to store information about regions used
#' for indicating regions for pileup
#'
#' @field .extptr an externalptr to c-level samtools regidx
#' @field gr GRanges containing regions
#' @field open logical indicating if the index is open
#' @rdname index_bed
.RegionIndex <- setRefClass("RegionIndex",
                        fields = list(
                          .extptr = "externalptr",
                          gr = "GRanges",
                          open = "logical"
                        )
)
