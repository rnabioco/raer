#' Annotate known SNP positions
#'
#' @description This function will annotate a GRanges or rowRanges of
#' a SummarizedExperiment with SNP positions from a SNP package
#'
#' @param obj GRanges or SummarizedExperiment  object
#' @param dbsnp SNPlocs package, see available packages from
#' [BSgenome::available.SNPs()]
#' @param chrom only operate on a specified chromosome
#' @param col_to_aggr column from SNPlocs package to add to
#' input. If multiple SNPs overlap these values will be concatenated
#' as comma separated values.
#' @param drop If TRUE, remove sites overlap SNPs
#' @param ... For the generic, further arguments to pass to specific methods.
#' Unused for now.
#'
#' @return Either a GRanges or SummarizedExperiment object with
#' a new column "snp" added with information from "col_to_aggr"
#'
#' @examples
#' if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
#'     gr <- GRanges(rep("22", 10),
#'         IRanges(
#'             seq(10510077,
#'                 10610077,
#'                 by = 1000
#'             )[1:10],
#'             width = 250
#'         ),
#'         strand = "+"
#'     )
#'     annot_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
#' }
#' @seealso [SNPlocs.Hsapiens.dbSNP144.GRCh38](https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh38.html)
#' @export
annot_snps <- function(obj, ...) {
    UseMethod("annot_snps", obj)
}
