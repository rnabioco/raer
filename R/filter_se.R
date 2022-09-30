
#' Remove multiallelic sites
#'
#' @description Sites with multiple variant bases will be removed from the
#' SummarizedExperiment. The rowData() will be updated to include a new column,
#' "Var", containing the variant allele detected at each site.
#'
#' @param se SummarizedExperiment
#'
#' @examples
#' example(create_se, echo = FALSE)
#' remove_multiallelic(se)
#'
#' @rdname filter_se
#'
#' @importFrom stringr str_count
#' @export
remove_multiallelic <- function(se) {
  is_not_multiallelic <- apply(assay(se, "Var"), 1, function(x) {
    x <- unique(x[x != "-"])
    if (length(x) == 0 | length(x) >= 2) {
      return(NA)
    }
    stringr::str_count(x, ",") == 0
  })
  se <- se[which(is_not_multiallelic), ]
  rowData(se)$Var <- apply(assay(se, "Var"), 1, function(x) unique(x[x != "-"]))
  se
}


#' Extract regions surrounding splice sites
#'
#' @description This function will return intervals containing splice sites and
#' adjacent regions.
#'
#' @param txdb A TxDb object
#' @param slop The number of bases upstream and downstream of splice site to extract
#' @return
#' GenomicRanges object containing positions of splice sites, with flanking bases.
#' @examples
#' if(require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
#'   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'   res <- get_splice_sites(txdb)
#'   res[1:5]
#' }
#' @importFrom GenomicFeatures intronsByTranscript
#' @export
get_splice_sites <- function(txdb, slop = 4) {
  if (!is(txdb, "TxDb")) {
    stop("txdb must be a TxDb object")
  }

  gr <- GenomicFeatures::intronsByTranscript(txdb)
  usub <- unlist(gr)

  int_start <- GRanges(seqnames(usub),
    IRanges(start(usub) - slop,
      start(usub) + slop - 1),
    strand = strand(usub))
  int_end <- GRanges(seqnames(usub),
    IRanges(end(usub) - slop - 1,
      end(usub) + slop),
    strand = strand(usub))
  int_pos <- c(int_start, int_end)
  sort(int_pos)
}

#' Remove sites near splice sites
#'
#' @description This function will remove editing sites found in regions proximal
#' to annotated splice junctions.
#'
#' @param se SummarizedExperiment with editing sites
#' @param txdb A TxDb object
#' @param splice_site_dist distance to splice site
#' @param ignore.strand if TRUE do not consider strand when comparing editing sites to
#' splice sites
#'
#' @rdname filter_se
#' @examples
#' if(require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
#'   library(SummarizedExperiment)
#'   nrows <- 5; ncols <- 6
#'   counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#'   rowRanges <- GRanges(rep("chr1", 5),
#'                     IRanges(c(12055, 12174, 12194, 12719, 12889), width=1),
#'                     strand=rep("+", 5))
#'   colData <- DataFrame(Treatment=rep(c("adar_wt", "adar_ko"), 3),
#'                        row.names=LETTERS[1:6])
#'   rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
#'                            rowRanges=rowRanges, colData=colData)
#'
#'   se <- remove_splice_variants(rse, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'   se
#' }
#' @importFrom GenomicFeatures intronsByTranscript
#' @export
remove_splice_variants <- function(se, txdb,
                                   splice_site_dist = 4,
                                   ignore.strand = FALSE) {
  ss <- get_splice_sites(txdb, splice_site_dist)
  x <- rowRanges(se)
  fo <- findOverlaps(x, ss, type = "any", ignore.strand = ignore.strand)
  to_keep <- setdiff(1:length(x), unique(queryHits(fo)))
  se[to_keep, ]
}

#' Remove clustered sequence variants
#'
#' @description Sequence variants of multiple allele types (e.g. A->G, C->T)
#' in nearby regions can be due to misalignments. This function will remove
#' variants if multiple allele types are present within a given distance
#' in genomic coordinate space or transcriptome coordinate space.
#'
#' @param se SummarizedExperiment containing editing sites
#' @param txdb A TxDb object
#' @param regions One of "transcript" and/or "genome", coordinate space to use
#' to examine distances between variants.
#' @param variant_dist distance in nucleotides for determining clustered variants
#'
#' @examples
#' if(require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
#'   library(SummarizedExperiment)
#'   nrows <- 5; ncols <- 6
#'   counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#'   rowRanges <- GRanges(rep("chr1", 5),
#'                     IRanges(c(12055, 12174, 12194, 12719, 12889), width=1),
#'                     strand=rep("+", 5))
#'   mcols(rowRanges)$Var <- c("AG", "AT", "AC", "TC", "GC")
#'   colData <- DataFrame(Treatment=rep(c("adar_wt", "adar_ko"), 3),
#'                        row.names=LETTERS[1:6])
#'   rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
#'                            rowRanges=rowRanges, colData=colData)
#'
#'   se <- remove_clustered_variants(rse, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'   se
#' }
#' @rdname filter_se
#' @return
#' SummarizedExperiment with sites removed from object dependent on filtering applied.
#'
#' @importFrom GenomicFeatures mapToTranscripts
#' @export
remove_clustered_variants <- function(se, txdb,
                                      regions = c("transcript", "genome"),
                                      variant_dist = 100) {

  if (!is(txdb, "TxDb")) {
    stop("txdb must be a TxDb object")
  }

  if (length(setdiff(regions, c("transcript", "genome"))) > 0) {
    stop("only transcript and/or genome are valid arguments for region")
  }

  x <- rowRanges(se)
  if ("genome" %in% regions) {
    fo <- findOverlaps(x, x + variant_dist)
    vars <- split(x[subjectHits(fo)]$Var, queryHits(fo))
    to_keep <- names(vars)[unlist(lapply(vars,
      function(x) {
        length(unique(x)) == 1
      }))]

    x <- x[as.integer(to_keep)]
  }

  if ("transcript" %in% regions) {
    tx_sites <- mapToTranscripts(x, txdb)
    tx_sites$Var <- x[tx_sites$xHits]$Var
    tx_sites <- sort(tx_sites)
    slop <- trim(suppressWarnings(tx_sites + variant_dist))

    fo <- findOverlaps(tx_sites, slop)
    vars <- split(tx_sites[subjectHits(fo)]$Var, queryHits(fo))
    to_drop <- names(vars)[unlist(lapply(vars,
      function(x) {
        length(unique(x)) > 1
      }))]
    tx_sites <- tx_sites[as.integer(to_drop)]
    x <- x[setdiff(1:length(x), unique(tx_sites$xHits))]
  }
  se[names(x), ]
}
