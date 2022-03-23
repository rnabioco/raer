#' Annotate snps
#'
#' @param obj GRanges or SummarizedExperiment  object
#' @param dbsnp SNPlocs package, see available packages from
#' [BSgenome::available.SNPs()]
#' @param chrom only operate on a specified chromosome
#' @param col_to_aggr column from SNPlocs package to add to
#' input. If multiple SNPs overlap these values will be concatenated
#' as comma separated values.
#' @param drop If TRUE, remove sites overlap SNPs
#'
#' @returns Either a GRanges or SummarizedExperiment object with
#' a new column "snp" added with information from "col_to_aggr"
#' @export
annot_snps <- function(obj, ...) {
  UseMethod("annot_snps", obj)
}

#' Index Bed
#' @export
setGeneric("indexBed",
           function(file, ...) standardGeneric("indexBed"))

