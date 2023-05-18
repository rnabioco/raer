#' Filter out multi-allelic sites
#'
#' @description Remove sites with multiple variant bases from a
#'   `SummarizedExperiment`. `rowData()` gains a new column, `ALT`, that
#'   contains the variant allele detected at each site.
#'
#' @param se `SummarizedExperiment::SummarizedExperiment`
#'
#' @examples
#' data(rse_adar_ifn)
#' filter_multiallelic(rse_adar_ifn)
#'
#' @returns `SummarizedExperiment::SummarizedExperiment` with multiallelic sites
#' removed.  A new column,`ALT` will be added to `rowData()` indicating the
#' single allele present at the site.
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
        !grepl(",", x)
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
#'     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'     res <- get_splice_sites(txdb)
#'     res[1:5]
#' }
#'
#' @importFrom GenomicFeatures intronsByTranscript
#'
#' @export
get_splice_sites <- function(txdb, slop = 4) {
    if (!is(txdb, "TxDb")) {
        cli::cli_abort("txdb must be a TxDb object")
    }

    int_gr <- GenomicFeatures::intronsByTranscript(txdb)
    int_gr <- unlist(int_gr)

    int_start <- GRanges(seqnames(int_gr),
        IRanges(
            start(int_gr) - slop,
            start(int_gr) + slop - 1
        ),
        strand = strand(int_gr)
    )
    int_end <- GRanges(seqnames(int_gr),
        IRanges(
            end(int_gr) - slop + 1,
            end(int_gr) + slop
        ),
        strand = strand(int_gr)
    )
    int_pos <- c(int_start, int_end)
    sort(int_pos)
}

#' Filter out sites near splice sites
#'
#' @description Remove editing sites found in regions proximal to annotated
#'   splice junctions.
#'
#' @param rse `SummarizedExperiment::SummarizedExperiment` with editing sites
#' @param txdb `GenomicFeatures::TxDb`
#' @param splice_site_dist distance to splice site
#' @param ignore.strand if `TRUE`, ignore strand when comparing editing sites to
#'   splice sites
#'
#' @family se-filters
#'
#' @examples
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'     data(rse_adar_ifn)
#'     gr <- GRanges(c(
#'         "DHFR:310-330:-",
#'         "DHFR:410-415:-",
#'         "SSR3:100-155:-",
#'         "SSR3:180-190:-"
#'     ))
#'     gr$source <- "raer"
#'     gr$type <- "exon"
#'     gr$source <- NA
#'     gr$phase <- NA_integer_
#'     gr$gene_id <- c(1, 1, 2, 2)
#'     gr$transcript_id <- rep(c("1.1", "2.1"), each = 2)
#'     txdb <- makeTxDbFromGRanges(gr)
#'     filter_splice_variants(rse_adar_ifn, txdb)
#' }
#'
#' @returns `SummarizedExperiment::SummarizedExperiment` with sites
#' adjacent to splice sites removed.
#'
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom GenomeInfoDb keepSeqlevels
#' @export
filter_splice_variants <- function(rse, txdb,
    splice_site_dist = 4,
    ignore.strand = FALSE) {
    n_in <- nrow(rse)

    spl_sites <- get_splice_sites(txdb, splice_site_dist)
    shared_seqs <- intersect(
        seqnames(seqinfo(rse)),
        seqnames(seqinfo(spl_sites))
    )
    if (length(shared_seqs) == 0) {
        cli::cli_abort("No shared seqnames found between txdb and rse")
    }
    spl_sites <- spl_sites[seqnames(spl_sites) %in% shared_seqs, ]
    spl_sites <- GenomeInfoDb::keepSeqlevels(spl_sites, shared_seqs)
    x <- rowRanges(rse)
    fo <- findOverlaps(x, spl_sites,
        type = "any",
        ignore.strand = ignore.strand
    )
    to_keep <- setdiff(seq_along(x), unique(queryHits(fo)))

    n_filt <- length(to_keep)
    cli::cli_alert_info(
        c(
            "{.fun filter_splice_variants}: removed {.val {n_in - n_filt}} sites",
            " from {.val {n_in}} ({.val {n_filt}} remain)"
        )
    )

    rse[to_keep, ]
}

#' Filter out clustered sequence variants
#'
#' @description Sequence variants of multiple allele types (e.g., `A>G`, `A>C`)
#'   in nearby regions can be due to mis-alignment. Remove variants if multiple
#'   allele types are present within a given distance in genomic or
#'   transcriptome coordinate space.
#'
#' @param rse `SummarizedExperiment::SummarizedExperiment` containing editing sites
#' @param txdb `GenomicFeatures::TxDb`
#' @param regions One of `transcript` or `genome`, specifying the coordinate
#'   system for calculating distances between variants.
#' @param variant_dist distance in nucleotides for determining clustered
#'   variants
#'
#' @examples
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'     data(rse_adar_ifn)
#'     gr <- GRanges(c(
#'         "DHFR:310-330:-",
#'         "DHFR:410-415:-",
#'         "SSR3:100-155:-",
#'         "SSR3:180-190:-"
#'     ))
#'     gr$source <- "raer"
#'     gr$type <- "exon"
#'     gr$source <- NA
#'     gr$phase <- NA_integer_
#'     gr$gene_id <- c(1, 1, 2, 2)
#'     gr$transcript_id <- rep(c("1.1", "2.1"), each = 2)
#'     txdb <- makeTxDbFromGRanges(gr)
#'
#'     rse <- filter_multiallelic(rse_adar_ifn)
#'
#'     filter_clustered_variants(rse, txdb)
#' }
#'
#' @family se-filters
#'
#' @return `SummarizedExperiment::SummarizedExperiment` with sites removed from
#'   object dependent on filtering applied.
#'
#' @importFrom GenomicFeatures mapToTranscripts
#' @export
filter_clustered_variants <- function(rse, txdb,
    regions = c("transcript", "genome"),
    variant_dist = 100) {
    if (!is(txdb, "TxDb")) {
        cli::cli_abort("txdb must be a TxDb object")
    }

    if (length(setdiff(regions, c("transcript", "genome"))) > 0) {
        cli::cli_abort("only transcript and/or genome are valid arguments for region")
    }

    n_in <- nrow(rse)

    x <- rowRanges(rse)

    if ("genome" %in% regions) {
        x_extend <- trim(suppressWarnings(x + variant_dist))
        fo <- findOverlaps(x, x_extend)
        fo_vars <- paste0(x[subjectHits(fo)]$REF, x[subjectHits(fo)]$ALT)
        vars <- split(fo_vars, queryHits(fo))
        to_keep <- names(vars)[unlist(lapply(
            vars,
            function(x) {
                length(unique(x)) == 1
            }
        ))]

        gn_keep <- as.integer(to_keep)
    } else {
        gn_keep <- seq_along(x)
    }

    if ("transcript" %in% regions) {
        x_tx <- x
        shared_seqs <- intersect(seqnames(x), seqnames(seqinfo(txdb)))
        if (length(shared_seqs) == 0) {
            cli::cli_abort("No shared seqnames found between txdb and rse")
        }
        x_tx <- x
        x_tx$id <- seq_along(x)
        x_tx <- x_tx[seqnames(x_tx) %in% shared_seqs]
        x_tx <- keepSeqlevels(x_tx, shared_seqs)
        tx_sites <- mapToTranscripts(x_tx,
            txdb,
            extractor.fun = GenomicFeatures::exonsBy
        )
        tx_sites$Var <- paste0(
            x_tx[tx_sites$xHits]$REF,
            x_tx[tx_sites$xHits]$ALT
        )
        tx_sites$id <- x_tx[tx_sites$xHits]$id
        tx_sites <- sort(tx_sites)
        tx_extend <- trim(suppressWarnings(tx_sites + variant_dist))

        fo <- findOverlaps(tx_sites, tx_extend)
        fo_vars <- tx_sites[subjectHits(fo)]$Var
        vars <- split(fo_vars, queryHits(fo))
        to_drop <- names(vars)[unlist(lapply(
            vars,
            function(x) {
                length(unique(x)) > 1
            }
        ))]
        tx_sites <- tx_sites[as.integer(to_drop)]
        tx_keep <- setdiff(seq_along(x), unique(tx_sites$id))
    } else {
        tx_keep <- seq_along(x)
    }

    x <- x[intersect(gn_keep, tx_keep), ]

    n_out <- length(x)

    cli::cli_alert_info(
        c(
            "{.fun filter_clustered_variants}: removed {.val {n_in - n_out}} sites",
            " from {.val {n_in}} ({.val {n_out}} remain)"
        )
    )

    rse[names(x), ]
}


#' Calculate confidence score for observing editing
#'
#' @description Calculate a confidence score based on a Bayesian inverse probability
#' model as described by Washburn et al. Cell Reports. 2015, and implemented
#' in the SAILOR pipeline.
#'
#' @param se `SummarizedExperiment::SummarizedExperiment` containing editing sites
#' @param edit_to edited base
#' @param edit_from non-edited base
#' @param per_sample if TRUE, calculate confidence per sample, otherwise edited
#' and non-edited counts will be summed across all samples.
#' @param exp_fraction Numeric, confidence margin parameter for
#'
#' @examples
#' data(rse_adar_ifn)
#' calc_confidence(rse_adar_ifn)
#' calc_confidence(rse_adar_ifn, per_sample = TRUE)
#'
#' @return `SummarizedExperiment::SummarizedExperiment` with either a new assay
#' or rowData column named "confidence" depending on whether confidence is
#'  calculated `per_sample`.
#'
#' @references
#' Washburn MC, Kakaradov B, Sundararaman B, Wheeler E, Hoon S, Yeo GW, Hundley HA. The dsRBP and inactive editor ADR-1 utilizes dsRNA binding to regulate A-to-I RNA editing across the C. elegans transcriptome. Cell Rep. 2014 Feb 27;6(4):599-607. doi: 10.1016/j.celrep.2014.01.011. Epub 2014 Feb 6. PMID: 24508457; PMCID: PMC3959997.
#'
#' SAILOR pipeline: https://github.com/YeoLab/sailor
#' @importFrom stats pbeta
#' @export
calc_confidence <- function(se,
    edit_to = "G",
    edit_from = "A",
    per_sample = FALSE,
    exp_fraction = 0.01) {
    if (length(exp_fraction) != 1) {
        cli::cli_abort("exp_fraction must be numeric(1)")
    }
    edit_to <- paste0("n", edit_to)
    edit_from <- paste0("n", edit_from)
    alt <- assay(se, edit_to)
    ref <- assay(se, edit_from)
    if (per_sample) {
        nc <- ncol(se)
        res <- vapply(seq_len(nc), function(i) {
            1 - pbeta(exp_fraction, alt[, i], ref[, i])
        }, FUN.VALUE = numeric(nrow(se)))
        colnames(res) <- colnames(se)
        assays(se)$confidence <- res
    } else {
        alt <- rowSums(alt)
        ref <- rowSums(ref)
        res <- 1 - pbeta(exp_fraction, alt, ref)
        rowData(se)$confidence <- res
    }
    se
}
