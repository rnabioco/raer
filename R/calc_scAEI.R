#' Calculate the Adenosine Editing Index (AEI) in single cells
#'
#' @description The Adenosine Editing Index describes the magnitude of A-to-I
#'   editing in a sample. The index is a weighted average of editing events (G
#'   bases) observed at A positions. The vast majority A-to-I editing occurs in
#'   ALU elements in the human genome, and these regions have a high A-to-I
#'   editing signal compared to other regions such as coding exons. This
#'   function will perform pileup at specified repeat regions and return a
#'   summary AEI metric.
#'
#' @references Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide
#' quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods
#' 16, 1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9
#'
#' @param bamfiles a path to a BAM file (for 10x libraries), or a vector of paths
#' to BAM files (smart-seq2). Can be supplied as a character vector, [BamFile], or
#' [BamFileList].
#' @param sites a GRanges object produced by [get_scAEI_sites] containing sites to process.
#' @param cell_barcodes A character vector of single cell barcodes to process. If
#' processing multiple BAM files (e.g. smart-seq-2), provide a character vector
#' of unique identifiers for each input BAM, to name each BAM file in the output files.
#' @param param object of class [FilterParam()] which specify various filters to
#'   apply to reads and sites during pileup.
#' @param edit_from This should correspond to the base
#'  (`A`, `C`, `G`, `T`) you expect in the reference. Ex. for A to I
#'   editing events, this would be `A`.
#' @param edit_to This should correspond to the base
#'  (`A`, `C`, `G`, `T`) you expect in an edited site. Ex. for A
#'   to I editing events, this would be `G`.
#' @param output_directory Output directory for `nRef` and `nAlt` sparseMatrix files.
#' If NULL, a temporary directory will be used.
#' @param return_sce if `TRUE`, data is returned as a SingleCellExperiment, if
#' `FALSE` a [DataFrame] containing computed AEI values will be returned.
#' @param ... additional arguments to [pileup_cells]
#'
#' @returns A `SingleCellExperiment` or a `DataFrame` containing computed AEI values,
#' count of editing events (n_alt), and count of referenc events (n_ref) per cell.
#'
#' @rdname calc_scAEI
#' @export
calc_scAEI <- function(bamfiles, sites, cell_barcodes, param = FilterParam(),
                       edit_from = "A", edit_to = "G",
                       output_dir = NULL, return_sce = FALSE,
                       ...) {

    if(is.null(output_dir)) {
        output_dir <- tempdir()
        on.exit(unlink(output_dir, recursive = TRUE))
    }

    # if unstranded, only query w.r.t + strand
    if(!param@library_type %in% c(1, 2)) {
        is_minus <- strand(sites) == "-"
        sites[is_minus]$REF <- comp_bases(edit_from)
        sites[is_minus]$ALT <- comp_bases(edit_to)
        strand(sites[is_minus]) <- "+"
    }

    aei_sce <- pileup_cells(
        bamfiles,
        sites,
        cell_barcodes,
        output_directory = output_dir,
        return_sce = TRUE,
        param = param,
        ...
    )

    n_alt <- Matrix::colSums(assay(aei_sce, "nAlt"))
    n_ref <- Matrix::colSums(assay(aei_sce, "nRef"))
    aei <- 100 * (n_alt / (n_alt + n_ref))
    res <- DataFrame(row.names = colnames(aei_sce),
                     AEI = aei,
                     n_alt = n_alt,
                     n_ref = n_ref)

    if(return_sce) {
        colData(aei_sce) <- res
        return(aei_sce)
    }
    res
}

#' Calculate the Adenosine Editing Index (AEI) in single cells
#'
#' @description The Adenosine Editing Index describes the magnitude of A-to-I
#'   editing in a sample. The index is a weighted average of editing events (G
#'   bases) observed at A positions. The vast majority A-to-I editing occurs in
#'   ALU elements in the human genome, and these regions have a high A-to-I
#'   editing signal compared to other regions such as coding exons. This
#'   function will perform pileup at specified repeat regions and return a
#'   summary AEI metric.
#'
#' @references Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide
#' quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods
#' 16, 1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9
#'
#' @param fasta_fn a path to a BAM file (for 10x libraries), or a vector of paths
#' to BAM files (smart-seq2). Can be supplied as a character vector, [BamFile], or
#' [BamFileList].
#' @param genes  A [GRanges] object with gene coordinates.Alternatively a [TxDb] object,
#' which if supplied, will be used extract gene coordinates.
#' @param alus [GRanges] with repeat regions to query for calculating the AEI,
#' typically ALU repeats. The strand of the supplied intervals will be ignored for
#' defining repeat regions.
#' @param edit_from This should correspond to the base
#'  (`A`, `C`, `G`, `T`) you expect in the reference. Ex. for A to I
#'   editing events, this would be `A`.
#' @param edit_to This should correspond to the base
#'  (`A`, `C`, `G`, `T`) you expect in an edited site. Ex. for A
#'   to I editing events, this would be `G`.
#'
#' @rdname calc_scAEI
#' @export
get_scAEI_sites <- function(fasta_fn, genes, alus,
                            edit_from = "A", edit_to = "G") {
    if(is(genes, "TxDb")) {
        genes <- suppressMessages(GenomicFeatures::genes(genes))
    }
    if(!is(genes, "GRanges")) {
        cli::cli_abort("genes must be a GRanges or TxDb object")
    }
    if(!is(alus, "GRanges")) {
        cli::cli_abort("alus must be a GRanges or TxDb object")
    }

    mcols(alus) <- NULL

    # store an integer id to allow user to track sites corresponding to supplied alu
    alus$id <- seq_along(alus)
    gn_alus <- prep_genic_alu_regions(genes, alus)
    aei_sites <- get_aei_site_positions(gn_alus,
                                        fasta_fn,
                                        edit_from)

    aei_sites$ALT <- edit_to

    # Rle encode to save memory
    aei_sites$ALT <- S4Vectors::Rle(aei_sites$ALT)
    aei_sites$id <- S4Vectors::Rle(aei_sites$id)
    aei_sites
}

prep_genic_alu_regions <- function(genes_gr, alu_gr) {
    # annotate gene region strandedness
    gns_ovl <- disjoin(genes_gr, with.revmap=TRUE, ignore.strand=TRUE)
    gn_strands <- unique(extractList(strand(genes_gr), gns_ovl$revmap))
    gn_strands[lengths(gn_strands) > 1] <- "*"
    strand(gns_ovl) <- unlist(gn_strands)

    # annotate genic alus strandedness based on gene
    hits <- findOverlaps(alu_gr, gns_ovl, ignore.strand = TRUE)
    genic_alus <- alu_gr[queryHits(hits), ]
    alu_gn_strands <- unique(strand(extractList(gns_ovl, hits)))
    alu_gn_strands[lengths(alu_gn_strands) > 1] <- "*"
    strand(genic_alus) <- unlist(alu_gn_strands[queryHits(hits), ])

    ss_alus <- genic_alus[strand(genic_alus) != "*"]
    ds_alus <- genic_alus[strand(genic_alus) == "*"]
    ss_alus$gene_strand <- "defined"

    ds_alus <- rep(ds_alus, each = 2)
    strand(ds_alus) <- c("+", "-")
    ds_alus$gene_strand <- "ambiguous"

    genic_alus <- unique(c(ss_alus, ds_alus))
    genic_alus$gene_strand <- S4Vectors::Rle(genic_alus$gene_strand)
    genic_alus
}

aei_site_positions <- function(base, seqs, gr) {
    base_pos <- Biostrings::vmatchPattern(base, seqs)
    tmp_gr <- gr
    alu_base_gr <- rep(tmp_gr, elementNROWS(base_pos))

    og_starts <- start(alu_base_gr)
    og_ends <- end(alu_base_gr)
    minus_as <- as.logical(strand(alu_base_gr) == "-")
    si <- unlist(Biostrings::startIndex(base_pos))
    start(alu_base_gr) <- og_starts + si - 1
    start(alu_base_gr)[minus_as] <- og_ends[minus_as] - si[minus_as] + 1
    end(alu_base_gr) <- start(alu_base_gr)
    alu_base_gr
}

get_aei_site_positions <- function(gr, fasta, base) {
    if(any(strand(gr) == "*")) {
        cli::cli_abort("strand must be set to + or -")
    }
    alu_seqs <- Rsamtools::scanFa(fasta, gr)
    minus_alus <- strand(gr) == "-"
    alu_seqs[minus_alus] <- Biostrings::reverseComplement(alu_seqs[minus_alus])
    res <- aei_site_positions(base, alu_seqs, gr)
    res$REF <- base
    res$REF <- S4Vectors::Rle(res$REF)
    res
}
