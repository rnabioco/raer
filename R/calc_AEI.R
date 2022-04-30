#' Calculate the Adenosine Editing Index (AEI)
#'
#' @description The Adenosine Editing Index describes the magnitude of A-to-I editing
#' in a sample. The index is a weighted average of editing events (G bases) observed
#' at A positions. The vast majority A-to-I editing occurs in ALU elements in the human
#' genome, and these regions have a high A-to-I editing signal compared to other regions
#' such as coding exons. This function will perform pileup at specified repeat regions and
#' return a summary AEI metric.
#'
#' @references
#' Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods 16, 1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9
#'
#'
#' @param bam_fn bam file
#' @param fasta_fn fasta
#' @param alu_ranges GRanges with regions to query for calculating the AEI,
#' typically ALU repeats.
#' @param txdb A txdb object, if supplied,  will be used to subset the alu_ranges to
#' those found overlapping genes.
#' @param snp_db a SNP package, if supplied, will be used to exclude polymorphic positions
#' prior to calculating the AEI.
#' @param BPPARAM A [BiocParallelParam] object for specifying parallel options for
#' operating over chromosomes.
#'
#' @returns A named list with the AEI index computed for all allelic combinations.
#' If correctly computed the signal from the A_G index should be higher than other
#' alleles (T_C), which are most likely derived from noise or polymorphisms.
#'
#' @examples
#' bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#' calc_AEI(bamfn, fafn)
#'
#' @importFrom BiocParallel bpstop bplapply SerialParam
#' @importFrom GenomicFeatures genes
#' @importFrom rtracklayer export
#' @importFrom Rsamtools scanBamHeader
#' @export
calc_AEI <- function(bam_fn,
                     fasta_fn,
                     alu_ranges = NULL,
                     txdb = NULL,
                     snp_db = NULL,
                     BPPARAM = SerialParam()){

  chroms <- names(Rsamtools::scanBamHeader(bam_fn)[[1]]$targets)

  if(!is.null(alu_ranges)){
    alu_bed_fn <- tempfile(fileext = ".bed")

    if(!is.null(txdb)){
    gene_gr <- GenomicFeatures::genes(txdb)
    alu_ranges <- GRanges::subsetByOverlaps(alu_ranges, gene_gr, ignore.strand = TRUE)
    alu_ranges <- GRanges::reduce(alu_ranges)
    }
    rtracklayer::export(alu_ranges, alu_bed_fn)
    chroms <- intersect(chroms, as.character(unique(seqnames(alu_ranges))))
  } else {
    alu_bed_fn <- NULL
  }

  aei <- bplapply(seq_along(chroms), function(i){
    start <- Sys.time()
    message("\tworking on: ", chroms[i], " time: ", Sys.time())
    plp <- get_pileup(bam_fn,
                      fafile = fasta_fn,
                      bedfile = alu_bed_fn,
                      chrom = chroms[i],
                      min_reads = 1,
                      min_base_qual = 30,
                      min_mapq = 255,
                      library_type = c("fr-first-strand"),
                      event_filters = c(5, 5, 0, 0, 0, 0, 0),
                      only_keep_variants = FALSE)
    message("\tcompleted in : ", Sys.time() - start)

    if(!is.null(snp_db)){
      plp <- annot_snps(plp, snp_db)
      plp <- plp[plp$RefSNP_id != ""]
    }

    bases <- c("A", "T", "C", "G")
    var_list <- list()
    for(i in seq_along(bases)){
      rb <- bases[i]
      other_b <- setdiff(bases, rb)
      j <- plp[plp$Ref == rb]
      for(k in seq_along(other_b)){
        ab <- other_b[k]
        id <- paste0(rb, "_", ab)
        n_alt <- sum(mcols(j)[[paste0("n", ab)]])
        n_ref <- sum(mcols(j)[[paste0("n", rb)]])
        var_list[[id]]  <- c(alt = n_alt,
                             ref = n_ref,
                             prop = 0)
      }
    }
    var_list

  }, BPPARAM = BPPARAM)

  bpstop(BPPARAM)
  unlink(alu_bed_fn)

  names(aei) <- chroms
  aei_res <- lapply(seq_along(aei), function(x) {
    vals <- aei[[x]]
    id <- names(aei)[x]
    xx <- as.data.frame(t(do.call(data.frame, vals)))
    xx$allele <- rownames(xx)
    xx$chrom <- id
    rownames(xx) <- NULL
    xx})

  aei_res <- do.call(rbind, aei_res)
  aei_res <- split(aei_res, aei_res$allele)
  aei_res <- lapply(aei_res, function(x) 100 * (sum(x$alt) / (sum(x$ref) + sum(x$alt))))
  aei_res
}
