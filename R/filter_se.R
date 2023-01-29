#' Filter out multi-allelic sites
#'
#' @description Remove sites with multiple variant bases from a
#'   `SummarizedExperiment`. `rowData()` gains a new column, `Var`, that
#'   contains the variant allele detected at each site.
#'
#' @param se `SummarizedExperiment::SummarizedExperiment`
#'
#' @examples
#' filter_multiallelic(rse_adar_ifn)
#'
#' @family se-filters
#'
#' @export
filter_multiallelic <- function(se) {
  is_not_multiallelic <- apply(assay(se, "Var"), 1, function(x) {
    x <- unique(x[x != "-"])
    if (length(x) == 0 | length(x) >= 2) {
      return(NA)
    }
    !grepl(',', x)
  })
  se <- se[which(is_not_multiallelic), ]
  rowData(se)$Var <- apply(assay(se, "Var"), 1, function(x) unique(x[x != "-"]))
  se
}

#' Extract regions surrounding splice sites
#'
#' @description Find intervals containing splice sites and their adjacent
#'   regions.
#'
#' @param txdb `GenomicFeatures::TxDb`
#' @param slop The number of bases upstream and downstream of splice site to
#'   extract
#' @return `GenomicRanges::GRanges` containing positions of splice sites, with
#' flanking bases.
#'
#' @examples
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'   res <- get_splice_sites(txdb)
#'   res[1:5]
#' }
#'
#' @importFrom GenomicFeatures intronsByTranscript
#'
#' @export
get_splice_sites <- function(txdb, slop = 4) {
  if (!is(txdb, "TxDb")) {
    stop("txdb must be a TxDb object")
  }

  gr <- GenomicFeatures::intronsByTranscript(txdb)
  usub <- unlist(gr)

  int_start <- GRanges(seqnames(usub),
    IRanges(
      start(usub) - slop,
      start(usub) + slop - 1
    ),
    strand = strand(usub)
  )
  int_end <- GRanges(seqnames(usub),
    IRanges(
      end(usub) - slop - 1,
      end(usub) + slop
    ),
    strand = strand(usub)
  )
  int_pos <- c(int_start, int_end)
  sort(int_pos)
}

#' Filter out sites near splice sites
#'
#' @description Remove editing sites found in regions proximal to annotated
#'   splice junctions.
#'
#' @param se `SummarizedExperiment::SummarizedExperiment` with editing sites
#' @param txdb `GenomicFeatures::TxDb`
#' @param splice_site_dist distance to splice site
#' @param ignore.strand if `TRUE`, ignore strand when comparing editing sites to
#'   splice sites
#'
#' @family se-filters
#'
#' @examples
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'   library(SummarizedExperiment)
#'   nrows <- 5
#'   ncols <- 6
#'   counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#'   rowRanges <- GRanges(rep("chr1", 5),
#'     IRanges(c(12055, 12174, 12194, 12719, 12889), width = 1),
#'     strand = rep("+", 5)
#'   )
#'   colData <- DataFrame(
#'     Treatment = rep(c("adar_wt", "adar_ko"), 3),
#'     row.names = LETTERS[1:6]
#'   )
#'   rse <- SummarizedExperiment(
#'     assays = SimpleList(counts = counts),
#'     rowRanges = rowRanges, colData = colData
#'   )
#'
#'   se <- filter_splice_variants(rse, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'   se
#' }
#'
#' @importFrom GenomicFeatures intronsByTranscript
#'
#' @export
filter_splice_variants <- function(se, txdb,
                                   splice_site_dist = 4,
                                   ignore.strand = FALSE) {
  ss <- get_splice_sites(txdb, splice_site_dist)
  x <- rowRanges(se)
  fo <- findOverlaps(x, ss, type = "any", ignore.strand = ignore.strand)
  to_keep <- setdiff(1:length(x), unique(queryHits(fo)))
  se[to_keep, ]
}

#' Filter out clustered sequence variants
#'
#' @description Sequence variants of multiple allele types (e.g., `A>G`, `A>C`)
#'   in nearby regions can be due to mis-alignment. Remove variants if multiple
#'   allele types are present within a given distance in genomic or
#'   transcriptome coordinate space.
#'
#' @param se `SummarizedExperiment::SummarizedExperiment` containing editing sites
#' @param txdb `GenomicFeatures::TxDb`
#' @param regions One of `transcript` or `genome`, specifying the coordinate
#'   system for calculating distances between variants.
#' @param variant_dist distance in nucleotides for determining clustered
#'   variants
#'
#' @examples
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'   library(SummarizedExperiment)
#'   nrows <- 5
#'   ncols <- 6
#'   counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#'   rowRanges <- GRanges(rep("chr1", 5),
#'     IRanges(c(12055, 12174, 12194, 12719, 12889), width = 1),
#'     strand = rep("+", 5)
#'   )
#'   mcols(rowRanges)$Var <- c("AG", "AT", "AC", "TC", "GC")
#'   colData <- DataFrame(
#'     Treatment = rep(c("adar_wt", "adar_ko"), 3),
#'     row.names = LETTERS[1:6]
#'   )
#'   rse <- SummarizedExperiment(
#'     assays = SimpleList(counts = counts),
#'     rowRanges = rowRanges, colData = colData
#'   )
#'
#'   se <- filter_clustered_variants(rse, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'   se
#' }
#'
#' @family se-filters
#'
#' @return `SummarizedExperiment::SummarizedExperiment` with sites removed from
#'   object dependent on filtering applied.
#'
#' @importFrom GenomicFeatures mapToTranscripts
#' @export
filter_clustered_variants <- function(se, txdb,
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
    to_keep <- names(vars)[unlist(lapply(
      vars,
      function(x) {
        length(unique(x)) == 1
      }
    ))]

    x <- x[as.integer(to_keep)]
  }

  if ("transcript" %in% regions) {
    tx_sites <- mapToTranscripts(x, txdb)
    tx_sites$Var <- x[tx_sites$xHits]$Var
    tx_sites <- sort(tx_sites)
    slop <- trim(suppressWarnings(tx_sites + variant_dist))

    fo <- findOverlaps(tx_sites, slop)
    vars <- split(tx_sites[subjectHits(fo)]$Var, queryHits(fo))
    to_drop <- names(vars)[unlist(lapply(
      vars,
      function(x) {
        length(unique(x)) > 1
      }
    ))]
    tx_sites <- tx_sites[as.integer(to_drop)]
    x <- x[setdiff(1:length(x), unique(tx_sites$xHits))]
  }
  se[names(x), ]
}

#' Find regions with oligodT mispriming
#'
#' @description OligodT will prime at A-rich regions in an RNA. Reverse transcription
#' from these internal priming sites will install an oligodT sequence at the 3' end
#' of the cDNA. Sequence variants within these internal priming sites are enriched
#' for variants converting the genomic sequence to the A encoded by the oligodT primer.
#' Trimming poly(A) from the 3' ends of reads reduces but does not eliminate these signals
#'
#' This function will identify regions that are enriched for mispriming events. Reads
#' that were trimmed to remove poly(A) (encoded in the pa tag by 10x genomics) are
#' identified. The aligned 3' positions of these reads are counted, and sites passing
#' thresholds (at least 2 reads) are retained as possible sites of mispriming. Be default
#' regions 5 bases upstream and 20 bases downstream of these putative mispriming sites
#' are returned.
#'
#' @importFrom Rsamtools ScanBamParam BamFile open.BamFile close.Bamfile
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom IRanges grouplengths
#' @importFrom S4Vectors aggregate
#' @examples
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#' fa_fn <- raer_example("mouse_tiny.fasta")
#' find_mispriming_sites(bam_fn, fa_fn)
#'
#' @export
find_mispriming_sites <- function(bamfile, fafile, pos_5p = 5, pos_3p = 20,
                                  min_reads = 2, tag = "pa", tag_values = 6:300,
                                  n_reads_per_chunk = 1e6, verbose = TRUE){
  tg_lst <- list(tag_values)
  names(tg_lst) <- tag
  sbp <- Rsamtools::ScanBamParam(tagFilter = tg_lst,tag = tag)
  bf <- Rsamtools::BamFile(bamfile, yieldSize = n_reads_per_chunk)
  open(bf)
  pa_pks <- GRanges()
  repeat {
    # should only return for reads with pa tag set
    galn <- readGAlignments(bf, param = sbp)

    if (length(galn) == 0) break
    gr <- as(galn, "GRanges")
    if(verbose) {
      s_ivl <- gr[1]
      e_ivl <- gr[length(gr)]
      message("working on ", s_ivl, " to ", e_ivl)
    }

    # count # of overlapping reads
    ans <- merge_pa_peaks(gr)
    pa_pks <- c(pa_pks, ans)
  }
  close(bf)
  # merge again, handle edge cases between yieldsizes
  ans <- reduce(pa_pks, with.revmap = TRUE)
  mcols(ans) <- S4Vectors::aggregate(pa_pks, mcols(ans)$revmap,
                                     mean_pal = mean(pa_pks$mean_pal),
                                     n_reads = sum(pa_pks$n_reads),
                                     drop = FALSE)

  # keep reads above threshold, slop, and merge adjacent misprimed regions
  ans <- ans[ans$n_reads >= min_reads]
  ans <- resize(ans, pos_3p + width(ans))
  ans <- resize(ans, pos_5p + width(ans), fix = "end")

  res <- reduce(ans, with.revmap = TRUE)
  mcols(res) <- S4Vectors::aggregate(ans, mcols(res)$revmap,
                                     n_reads = sum(ans$n_reads),
                                     drop = FALSE)
  res$n_regions <- IRanges::grouplengths(res$grouping)
  res$grouping <- NULL
  res <- pa_seq_context(res, fafile)
  res
}

merge_pa_peaks <- function(gr) {
  # get 3' end of read
  start(gr[strand(gr) == "+"]) <- end(gr[strand(gr) == "+"])
  end(gr[strand(gr) == "-"]) <- start(gr[strand(gr) == "-"])

  # merge and count reads within merged ivls
  ans <- reduce(gr, with.revmap = TRUE)
  mcols(ans) <- S4Vectors::aggregate(gr, mcols(ans)$revmap,
                                     mean_pal = mean(gr$pa),
                                     drop = FALSE)
  mcols(ans)$n_reads <- IRanges::grouplengths(ans$grouping)
  ans
}

#' @importFrom Rsamtools FaFile scanFa
#' @importFrom Biostrings letterFrequency reverseComplement
pa_seq_context <- function(gr, fafile){
  fa <- Rsamtools::FaFile(fafile)
  seqs <- Rsamtools::scanFa(fa, gr)
  seqs[strand(gr) == "-"] <- Biostrings::reverseComplement(seqs[strand(gr) == "-"])
  a_prop <- Biostrings::letterFrequency(seqs, "A") / width(gr)
  mcols(gr)$A_freq <- a_prop[, 1]
  gr
}
