#' Calculate the Adenosine Editing Index (AEI)
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
#' @param bamfiles character vector of paths to indexed bam files. If a named
#' character vector is supplied the names will be used in the output.
#' @param fasta fasta filename
#' @param alu_ranges [GRanges] with regions to query for
#'   calculating the AEI, typically ALU repeats.
#' @param txdb A [TxDb] object, if supplied, will be used to subset the
#'   alu_ranges to those found overlapping genes. Alternatively a [GRanges]
#'   object with gene coordinates.  If the `library_type`, specified by
#'   `FilterParam`, is `unstranded` then the [TxDb] will
#'   be used to correct the strandness relative to the reference and is a
#'   required
#'   parameter.
#' @param snp_db either a [SNPlocs], [GPos], or [GRanges] object. If supplied,
#'   will be used to exclude polymorphic positions prior to calculating the AEI.
#'   If `calc_AEI()` will be used many times, one will save time by first
#'   identifying SNPs that overlap the supplied `alu_ranges`, and passing these
#'   as a [GRanges] to `snp_db` rather than supplying all known SNPs (see
#'   [get_overlapping_snps()]).
#' @param param object of class [FilterParam()] which specify various
#'   filters to apply to reads and sites during pileup.
#' @param BPPARAM A [BiocParallelParam] object for specifying parallel options
#'   for operating over chromosomes.
#' @param verbose report progress on each chromosome?
#'
#' @returns A named list containing:
#'   - `AEI`: a matrix of AEI index values computed for all allelic
#'   combinations, one row for each supplied bam file.
#'   - `AEI_per_chrom`: a data.frame containing values computed for each
#'   chromosome
#'
#' @examples
#' suppressPackageStartupMessages(library(Rsamtools))
#'
#' bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
#' bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
#' bams <- c(bamfn, bam2fn)
#' names(bams) <- c("ADAR1KO", "WT")
#'
#' fafn <- raer_example("human.fasta")
#' mock_alu_ranges <- scanFaIndex(fafn)
#'
#' res <- calc_AEI(bams, fafn, mock_alu_ranges)
#' res$AEI
#'
#' @importFrom BiocParallel bpstop bpmapply SerialParam
#' @importFrom GenomicFeatures genes
#' @importFrom rtracklayer export
#' @importFrom Rsamtools scanBamHeader
#' @importFrom IRanges subsetByOverlaps
#' @import S4Vectors
#' @import GenomicRanges
#'
#' @export
calc_AEI <- function(bamfiles,
    fasta,
    alu_ranges = NULL,
    txdb = NULL,
    snp_db = NULL,
    param = FilterParam(),
    BPPARAM = SerialParam(),
    verbose = FALSE) {
    if (!is(bamfiles, "BamFileList")) {
        if (is.null(names(bamfiles))) {
            names(bamfiles) <- bamfiles
        }
        bamfiles <- BamFileList(bamfiles)
    }

    chroms <- seqnames(seqinfo_from_header(bamfiles[[1]]))

    if (is.null(alu_ranges)) {
        cli::cli_alert_warning(c(
            "Querying the whole genome will be very ",
            "memory intensive and inaccurate.\n",
            "Consider supplying a GRanges object with ALU",
            "or related repeats for your species "
        ))
    }

    genes_gr <- NULL
    alu_bed_fn <- NULL
    if (param@library_type == 0) {
        if (is.null(txdb)) {
            cli::cli_abort("txdb required for processing unstranded data")
        }

        if (is(txdb, "TxDb")) {
            genes_gr <- GenomicFeatures::genes(txdb)
        } else {
            genes_gr <- txdb
        }
    }

    if (!is.null(alu_ranges)) {
        if (is(alu_ranges, "GRanges")) {
            if (!is.null(txdb)) {
                if (is(txdb, "TxDb")) {
                    genes_gr <- GenomicFeatures::genes(txdb)
                } else {
                    genes_gr <- txdb
                }
                alu_ranges <- subsetByOverlaps(alu_ranges,
                    genes_gr,
                    ignore.strand = TRUE
                )
                alu_ranges <- reduce(alu_ranges)
            }
            chroms <- intersect(
                chroms,
                as.character(unique(seqnames(alu_ranges)))
            )

            alu_ranges <- alu_ranges[seqnames(alu_ranges) %in% chroms]
        } else {
            cli::cli_abort("unrecognized format for alu_ranges")
        }
    }

    snps <- NULL
    if (!is.null(snp_db)) {
        if (is(snp_db, "GRanges") || is(snp_db, "GPos")) {
            if (is(alu_ranges, "GRanges")) {
                snps <- subsetByOverlaps(snp_db, alu_ranges)
            } else {
                snps <- snp_db
            }
        } else if (is(snp_db, "ODLT_SNPlocs")) {
            if (is(alu_ranges, "GRanges")) {
                alu_style <- seqlevelsStyle(alu_ranges)
                snp_style <- seqlevelsStyle(snp_db)
                if (!any(alu_style %in% snp_style)) {
                    cli::cli_alert_warning(c(
                        "seqlevels style in supplied snps ({snp_style}) ",
                        "differs from alu_ranges ({alu_style}) ",
                        "attempting to coerce"
                    ))
                    seqlevelsStyle(alu_ranges) <- snp_style
                    snps <- snpsByOverlaps(snp_db, alu_ranges)
                    seqlevelsStyle(snps) <- alu_style
                } else {
                    snps <- snpsByOverlaps(snp_db, alu_ranges)
                }
            } else {
                cli::cli_abort(
                    "removing snps using a SNPloc package requires alu_ranges"
                )
            }
        } else {
            cli::cli_abort("unknown snpdb object type")
        }
        snp_lst <- split(snps, seqnames(snps))
        missing_chroms <- setdiff(chroms, names(snp_lst))
        if (length(missing_chroms) > 0) {
            missing_snp_lst <- lapply(missing_chroms, function(x) GRanges())
            names(missing_snp_lst) <- missing_chroms
            snp_lst <- c(snp_lst, missing_snp_lst)
        }

        snps <- snp_lst[chroms]
    }
    res <- list()
    for (i in seq_along(bamfiles)) {
        bam_fn <- bamfiles[i]
        aei <- .AEI_per_bam(
            bam_fn = bam_fn, fasta = fasta, chroms = chroms,
            alu_ranges = alu_ranges, snps = snps, param = param,
            genes_gr = genes_gr, verbose = verbose, BPPARAM = BPPARAM
        )
        res[[path(bam_fn)]] <- aei
    }

    aei <- do.call(rbind, lapply(res, "[[", 1))
    rownames(aei) <- names(bamfiles)

    aei_per_chrom <- lapply(seq_along(res), function(i) {
        x <- res[[i]][[2]]
        data.frame(x, row.names = NULL)
    })
    aei_per_chrom <- do.call(rbind, aei_per_chrom)

    list(AEI = aei, AEI_per_chrom = aei_per_chrom)
}


.AEI_per_bam <- function(
        bam_fn, fasta, chroms, alu_ranges, snps, param,
        genes_gr, verbose,
        BPPARAM = BiocParallel::SerialParam()) {
    if (is.null(snps)) {
        aei <- bpmapply(.calc_AEI_per_chrom,
            chroms,
            MoreArgs = list(
                bam_fn = bam_fn,
                fasta = fasta,
                alu_sites = alu_ranges,
                param = param,
                snp_gr = NULL,
                genes_gr = genes_gr,
                verbose = verbose
            ),
            BPPARAM = BPPARAM,
            SIMPLIFY = FALSE
        )
    } else {
        if (length(chroms) != length(snps)) {
            cli::cli_abort("issue subsetting snpDb and chromosomes")
        }
        aei <- bpmapply(.calc_AEI_per_chrom,
            chroms,
            snps,
            MoreArgs = list(
                bam_fn = bam_fn,
                fasta = fasta,
                alu_sites = alu_ranges,
                param = param,
                genes_gr = genes_gr,
                verbose = verbose
            ),
            BPPARAM = BPPARAM,
            SIMPLIFY = FALSE
        )
    }
    bpstop(BPPARAM)

    names(aei) <- chroms
    aei <- lapply(seq_along(aei), function(x) {
        vals <- aei[[x]]
        id <- names(aei)[x]
        xx <- as.data.frame(t(do.call(data.frame, vals)))
        xx <- cbind(
            chrom = id,
            allele = rownames(xx),
            xx
        )
        rownames(xx) <- NULL
        xx
    })

    aei <- do.call(rbind, aei)
    aei <- split(aei, aei$allele)

    aei_summary <- lapply(aei, function(x) {
        100 * (sum(x$alt) / (sum(x$ref) + sum(x$alt)))
    })
    aei_summary <- t(do.call(rbind, aei_summary))

    aei_per_chrom <- do.call(rbind, aei)

    if (is.null(names(bam_fn))) {
        rownames(aei_summary) <- bam_fn
        aei_per_chrom$bam_file <- bam_fn
    } else {
        rownames(aei_summary) <- names(bam_fn)
        aei_per_chrom$bam_file <- names(bam_fn)
    }

    list(AEI = aei_summary, AEI_per_chrom = aei_per_chrom)
}

.calc_AEI_per_chrom <- function(bam_fn, fasta, alu_sites, chrom, param,
    snp_gr, genes_gr, verbose) {
    if (verbose) {
        start <- Sys.time()
        cli::cli_progress_step("working on: {chrom}")
    }
    param@min_depth <- 1L
    param@only_keep_variants <- FALSE

    plp <- pileup_sites(bam_fn,
        fasta = fasta,
        sites = alu_sites,
        chroms = chrom,
        param = param
    )

    if (length(snp_gr) > 0 && length(plp) > 0) {
        plp <- keepSeqlevels(plp, chrom)
        snp_gr <- keepSeqlevels(snp_gr, chrom)

        plp <- subsetByOverlaps(plp, snp_gr,
            invert = TRUE,
            ignore.strand = TRUE
        )
    }

    if (param@library_type == 0L) {
        plp <- correct_strand(plp, genes_gr)
    }

    bases <- c("A", "T", "C", "G")
    var_list <- list()
    for (i in seq_along(bases)) {
        rb <- bases[i]
        other_b <- setdiff(bases, rb)
        j <- plp[rowData(plp)$REF == rb]
        for (k in seq_along(other_b)) {
            ab <- other_b[k]
            id <- paste0(rb, "_", ab)
            n_alt <- sum(assays(j)[[paste0("n", ab)]][, 1])
            n_ref <- sum(assays(j)[[paste0("n", rb)]][, 1])
            var_list[[id]] <- c(
                alt = n_alt,
                ref = n_ref,
                prop = n_alt / (n_alt + n_ref)
            )
        }
    }
    var_list
}


#' Retrieve SNPs overlapping intervals
#'
#' @description This function will find SNPs overlapping supplied intervals
#'   using a SNPlocs package. The SNPs can be returned in memory (as GPos
#'   objects) or written to disk as a bed-file (optionally compressed).
#'
#' @param gr Intervals to query
#' @param snpDb A reference ot a SNPlocs database
#' @param output_file A path to an output file. If supplied the file can be
#'   optionally compressed by including a ".gz" suffix. If not supplied, SNPS
#'   will be returned as a [GenomicRanges::GPos] object
#'
#' @return GPos object containing SNPs overlapping supplied genomic intervals
#' @examples
#' if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
#'     gr <- GRanges(rep("22", 10),
#'         IRanges(seq(10510077, 10610077, by = 1000)[1:10], width = 250),
#'         strand = "+"
#'     )
#'     get_overlapping_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
#' }
#' @importFrom rtracklayer export
#' @importFrom BSgenome snpsByOverlaps
#'
#' @export
get_overlapping_snps <- function(gr,
    snpDb,
    output_file = NULL) {
    gr <- gr[seqnames(gr) %in% seqnames(snpDb)]

    # iterate through each contig, drop mcols (snpID) to reduce memory
    alu_snps <- vector("list", length = length(seqnames(snpDb)))
    for (i in seq_along(seqnames(snpDb))) {
        x <- seqnames(snpDb)[i]
        tmp_gr <- gr[seqnames(gr) == x]

        xx <- BSgenome::snpsByOverlaps(snpDb, tmp_gr)
        mcols(xx) <- NULL

        if (!is.null(output_file)) {
            rtracklayer::export(xx,
                output_file,
                append = TRUE,
                ignore.strand = TRUE
            )
        } else {
            alu_snps[[i]] <- xx
        }
    }
    if (!is.null(output_file)) {
        return(output_file)
    }
    unlist(as(alu_snps, "GRangesList"))
}

#' Apply strand correction using gene annotations
#'
#' @description Gene annotations are used to infer the likely strand of editing
#'   sites. This function will operate on unstranded datasets which have been
#'   processed using "unstranded" library type which reports variants
#'   with respect to the + strand for all sites. The strand of the editing site
#'   will be assigned the strand of overlapping features in the `genes_gr`
#'   object. Sites with no-overlap, or overlapping features with conflicting
#'   strands (+ and -) will be removed.
#'
#' @param rse RangedSummarizedExperiment object containing editing sites
#' processed with "unstranded" setting
#' @param genes_gr GRanges object containing reference features to annotate the
#'   strand of the editing sites.
#'
#' @return RangedSummarizedExperiment object containing pileup assays,
#' with strand corrected based on supplied genomic intervals.
#'
#' @examples
#' suppressPackageStartupMessages(library("GenomicRanges"))
#'
#' bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
#' fafn <- raer_example("human.fasta")
#' fp <- FilterParam(library_type = "unstranded")
#' rse <- pileup_sites(bamfn, fafn, param = fp)
#'
#' genes <- GRanges(c(
#'     "DHFR:200-400:+",
#'     "SPCS3:100-200:-",
#'     "SSR3:3-10:-",
#'     "SSR3:6-12:+"
#' ))
#'
#' correct_strand(rse, genes)
#'
#' @export
correct_strand <- function(rse, genes_gr) {
    if (length(rse) == 0) {
        return(rse)
    }

    stopifnot(all(strand(rse) == "+"))
    stopifnot("REF" %in% colnames(rowData(rse)))
    stopifnot("ALT" %in% names(assays(rse)))

    genes_gr$gene_strand <- strand(genes_gr)
    rse <- annot_from_gr(rse, genes_gr, "gene_strand", ignore.strand = TRUE)

    # drop non-genic and multi-strand (overlapping annotations)
    rse <- rse[!is.na(rowData(rse)$gene_strand), ]

    gs <- decode(rowData(rse)$gene_strand)
    n_strands <- lengths(regmatches(gs, gregexpr(",", gs)))
    rse <- rse[n_strands == 0, ]

    flip_rows <- as.vector(strand(rse) != rowData(rse)$gene_strand)

    rowData(rse)$REF[flip_rows] <- comp_bases(rowData(rse)$REF[flip_rows])

    if ("ALT" %in% colnames(rowData(rse))) {
        rowData(rse)$ALT[flip_rows] <- comp_bases(rowData(rse)$ALT[flip_rows])
    }

    # complement the ALT variants
    alts <- assay(rse, "ALT")[flip_rows, , drop = FALSE]
    comp_alts <- complement_variant_matrix(alts)
    assay(rse, "ALT")[flip_rows, ] <- comp_alts

    # complement the nucleotide counts by reordering the assays
    assays_to_swap <- c("nA", "nT", "nC", "nG")
    og_order <- rownames(rse)
    sites_to_swap <- rse[flip_rows, , drop = FALSE]

    to_swap <- assays(sites_to_swap)[assays_to_swap]
    to_swap <- to_swap[c(2, 1, 4, 3)]
    names(to_swap) <- assays_to_swap
    assays(sites_to_swap)[assays_to_swap] <- to_swap

    rse <- rbind(rse[!flip_rows, ], sites_to_swap)
    rse <- rse[og_order, ]

    strand(rse) <- rowData(rse)$gene_strand
    rowData(rse)$gene_strand <- NULL
    rse
}


# complement the ALT assay matrix,
# including multiallelics (e.g. c("A", "A,T", "C"))
complement_variant_matrix <- function(alt_mat) {
    stopifnot(all(!is.na(alt_mat)))

    # find multiallelics
    n_alleles <- lengths(regmatches(alt_mat, gregexpr(",", alt_mat)))

    # convert to a vector for processing
    alts <- as.vector(alt_mat)
    comp_alts <- alts

    # single allele sites
    comp_alts[n_alleles == 0] <- comp_bases(alts[n_alleles == 0])

    # split and complement multiallelics
    multialts <- alts[n_alleles > 0]
    multialts <- vapply(
        strsplit(multialts, ","), function(y) {
            paste0(comp_bases(y), collapse = ",")
        },
        character(1)
    )
    comp_alts[n_alleles > 0] <- multialts

    # get back to a matrix
    comp_alts <- matrix(comp_alts,
        nrow = nrow(alt_mat),
        ncol = ncol(alt_mat),
        dimnames = dimnames(alt_mat)
    )
    comp_alts
}
