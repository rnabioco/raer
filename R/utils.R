

#' Bin genomic regions
#' @description Genomic regions will be split first by chomosome if n < # of chromosomes,
#' otherwise chromosomes will be halved
#' @param regions named list with each entry from a chromosome and containing start and end coordinates
#' @param n desired number of splits
#
# bin_genome <- function(regions, n = ){
#
# }
#
# library(Rsamtools)
# library(magrittr)
# fa_regions <- scanFaIndex(fafile)
# seqlevels(fa_regions) <- seqlevels(fa_regions)[order(ranges(fa_regions), decreasing = TRUE)]
# fa_regions <- fa_regions[order(ranges(fa_regions), decreasing = TRUE)]
# fa_regions <- split(ranges(fa_regions), seqnames(fa_regions))
