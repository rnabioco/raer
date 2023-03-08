#' Sample data for RNA editing analysis
#'
#' A subset of data from an RNA-seq experiment to measure
#' the effects of IFN treatment of cell lines with wild-type
#' or ADAR1-KO
#'
#' @format ## `rse_adar_ifn`
#'
#' ```
#' class: RangedSummarizedExperiment
#' dim: 74 2
#' metadata(0):
#' assays(7): Var nRef ... nC nG
#' rownames(74): SSR3_102_- SSR3_125_- ... DHFR_430_- DHFR_513_-
#' rowData names(1): Ref
#' colnames(2): wt adar1_ko
#' colData names(1): sample
#' ```
#' @returns RangedSummarizedExperiment populated with pileup data
#' @source <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386593>
#' @references <https://pubmed.ncbi.nlm.nih.gov/29395325/>
#' @usage data(rse_adar_ifn)
"rse_adar_ifn"
