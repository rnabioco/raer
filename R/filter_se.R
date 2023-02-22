#' Filter out multi-allelic sites
#'
#' @description Remove sites with multiple variant bases from a
#'   `SummarizedExperiment`. `rowData()` gains a new column, `ALT`, that
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
  n_in <- nrow(se)
  is_not_multiallelic <- apply(assay(se, "ALT"), 1, function(x) {
    x <- unique(x[x != "-"])
    if (length(x) == 0 | length(x) >= 2) {
      return(NA)
    }
    !grepl(',', x)
  })
  se <- se[which(is_not_multiallelic), ]
  rowData(se)$ALT <- apply(assay(se, "ALT"), 1, function(x) unique(x[x != "-"]))

  n_filt <- sum(c(is.na(is_not_multiallelic), !is_not_multiallelic), na.rm = TRUE)
  cli::cli_alert_info(
    c(
      "{.fun filter_multiallelic}: removed {.val {n_filt}} sites",
      " from {.val {n_in}} ({.val {nrow(se)}} remain)"
    )
  )

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
#'   filter_splice_variants(
#'     rse_adar_ifn,
#'     TxDb.Hsapiens.UCSC.hg38.knownGene
#'   )
#' }
#'
#' @importFrom GenomicFeatures intronsByTranscript
#'
#' @export
filter_splice_variants <- function(se, txdb,
                                   splice_site_dist = 4,
                                   ignore.strand = FALSE) {
  n_in <- nrow(se)
  ss <- get_splice_sites(txdb, splice_site_dist)
  x <- rowRanges(se)
  fo <- findOverlaps(x, ss, type = "any", ignore.strand = ignore.strand)
  to_keep <- setdiff(1:length(x), unique(queryHits(fo)))

  n_filt <- length(to_keep)
  cli::cli_alert_info(
    c(
      "{.fun filter_splice_variants}: removed {.val {n_in - n_filt}} sites",
      " from {.val {n_in}} ({.val {n_filt}} remain)"
    )
  )

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
#'   rse <- filter_multiallelic(rse_adar_ifn)
#'
#'   filter_clustered_variants(
#'     rse,
#'     TxDb.Hsapiens.UCSC.hg38.knownGene
#'   )
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

  n_in <- nrow(se)

  x <- rowRanges(se)

  if ("genome" %in% regions) {
    fo <- findOverlaps(x, x + variant_dist)
    vars <- split(x[subjectHits(fo)]$ALT, queryHits(fo))
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
    tx_sites$ALT <- x[tx_sites$xHits]$ALT
    tx_sites <- sort(tx_sites)
    slop <- trim(suppressWarnings(tx_sites + variant_dist))

    fo <- findOverlaps(tx_sites, slop)
    vars <- split(tx_sites[subjectHits(fo)]$ALT, queryHits(fo))
    to_drop <- names(vars)[unlist(lapply(
      vars,
      function(x) {
        length(unique(x)) > 1
      }
    ))]
    tx_sites <- tx_sites[as.integer(to_drop)]
    x <- x[setdiff(1:length(x), unique(tx_sites$xHits))]
  }

  n_out <- length(x)

  cli::cli_alert_info(
    c(
      "{.fun filter_clustered_variants}: removed {.val {n_in - n_out}} sites",
      " from {.val {n_in}} ({.val {n_out}} remain)"
    )
  )

  se[names(x), ]
}
