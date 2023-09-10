#' Annotate known SNP positions
#'
#' @description This function will annotate a [GRanges] or the rowRanges of
#' a [SummarizedExperiment] with SNPs from a SNP package.
#'
#' @param obj GRanges or SummarizedExperiment  object
#' @param dbsnp SNPlocs package, see available packages from
#' [BSgenome::available.SNPs()]
#' @param chrom only operate on a specified chromosome
#' @param col_to_aggr column from SNPlocs package to add to
#' input. If multiple SNPs overlap these values will be concatenated
#' as comma separated values.
#' @param genome A BSgenome object, which if supplied, will be used to provide
#' additional `snp_ref_allele` and `snp_alt_alleles` columns containing the
#' reference and alt allele sequences, with respect to the positive strand.
#' Additionally the snp sequences will be checked against the allele at the site
#' if a column named `ALT` is present in object. The strand of the site will be
#' used to determine if the `ALT` allele needs to be complemented prior to
#' comparing against the SNP db (which always returns sequences w.r.t the
#' plus strand).
#' @param drop If TRUE, remove sites overlapping SNPs
#' @param RLE If TRUE, columns added will returned as [S4Vectors::Rle()] vectors
#' to reduce memory
#'
#' @param ... For the generic, further arguments to pass to specific methods.
#' Unused for now.
#'
#' @return Either a GRanges or SummarizedExperiment object with
#' a new column added with information from `col_to_aggr` and optionally
#' `snp_ref_allele`, `snp_alt_alleles`, and `snp_matches_site` annotations.
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
#'     genome(gr) <- "GRCh38.p2"
#'     annot_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
#' }
#' @seealso [SNPlocs.Hsapiens.dbSNP144.GRCh38](https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh38.html)
#' @export
annot_snps <- function(obj, ...) {
    UseMethod("annot_snps", obj)
}


#' @rdname annot_snps
#' @importFrom BSgenome snpsByOverlaps
#' @export
annot_snps.GRanges <- function(obj,
    dbsnp,
    chrom = NULL,
    col_to_aggr = "RefSNP_id",
    drop = FALSE,
    genome = NULL,
    RLE = TRUE,
    ...) {
    if (!is(dbsnp, "ODLT_SNPlocs")) {
        cli::cli_abort("supplied dbSNP not valid SNP package, please install SNP db")
    }

    if (any(col_to_aggr %in% colnames(mcols(obj)))) {
        cli::cli_abort("supplied column name in col_to_aggr already exists in input")
    }

    if (!is.null(chrom)) {
        sites <- obj[seqnames(obj) == chrom]
    } else {
        sites <- obj
    }

    in_style <- seqlevelsStyle(sites)
    snp_style <- seqlevelsStyle(dbsnp)
    if (!any(in_style %in% snp_style)) {
        cli::cli_alert_warning(c(
            "seqlevels style in supplied snps ({snp_style}) ",
            "differs from sites ({in_style}) ",
            "attempting to coerce"
        ))
    }

    seqlevelsStyle(sites) <- snp_style

    # returns GPos for each snp
    snps <- snpsByOverlaps(dbsnp, sites, genome = genome)

    # now annot if site is a SNP
    snp_overlaps <- findOverlaps(sites, snps, ignore.strand = TRUE)

    if (length(snp_overlaps) == 0) {
        mcols(sites)[col_to_aggr] <- NA
        seqlevelsStyle(sites) <- in_style
        return(sites)
    }

    if (is.null(genome)) {
        mcols(sites)[col_to_aggr] <- aggregate(
            snps,
            snp_overlaps,
            snp = unstrsplit(eval(parse(text = col_to_aggr)), ","),
            drop = FALSE
        )$snp

        if(RLE) {
            mcols(sites)[col_to_aggr] <- S4Vectors::Rle(mcols(sites)[[col_to_aggr]])
        }

    } else {
        cols_exist <- any(c("snp_ref_allele", "snp_alt_alleles") %in%
            colnames(mcols(obj)))
        if (cols_exist) {
            cli::cli_abort(
                c("snp_ref_allele or snp_alt_alleles columns already exist in input")
            )
        }

        # prevent no visible binding for global variable note
        alt_alleles <- ref_allele <- NULL

        snps$alt_alleles <- vapply(snps$alt_alleles,
                                   function(x) paste0(unique(x), collapse = ","),
                                   FUN.VALUE = character(1))

        snp_info <- aggregate(
            snps,
            snp_overlaps,
            snp = unstrsplit(eval(parse(text = col_to_aggr)), ","),
            ref_allele = unstrsplit(ref_allele),
            alt_alleles = unstrsplit(alt_alleles),
            drop = FALSE
        )

        snp_info$grouping <- NULL
        snp_cols <- c(col_to_aggr,
                      "snp_ref_allele",
                      "snp_alt_alleles")
        colnames(snp_info) <- snp_cols

        mcols(sites) <- cbind(mcols(sites), snp_info)

        if("ALT" %in% colnames(mcols(sites))) {
            snp_cols <- c(snp_cols, "snp_matches_site")
            mcols(sites)$snp_matches_site <- check_snp_match(sites)
        }

        if(RLE) {
            mcols(sites)[snp_cols] <- lapply(mcols(sites)[snp_cols],
                                             S4Vectors::Rle)
        }
    }

    seqlevelsStyle(sites) <- in_style

    if (drop) {
        sites <- sites[!is.na(sites$snp)]
    }

    sites
}


#' @rdname annot_snps
#' @export
annot_snps.SummarizedExperiment <- function(obj, ...) {
    gr <- rowRanges(obj)
    res <- annot_snps.GRanges(gr, ...)
    mcols(rowRanges(obj)) <- mcols(res)
    obj
}

#' Annotate sites using GRanges object
#'
#' @description Utility function to map annotations from GRanges to rowData of
#'   SummarizedExperiment or to mcols of GRanges object. If multiple features overlap then
#'   they will be concatenated as comma separated values.
#'
#' @param obj RangedSummarizedExperiment or GRanges object
#' @param gr GRanges with annotations to map to obj
#' @param cols_to_map character vector of columns from gr to map to obj. If the
#'   vector has names, the names will be the column names in the output obj
#' @param RLE If TRUE, columns added will returned as [S4Vectors::Rle()] vectors
#' to reduce memory
#' @param ... additional arguments to pass to [GenomicRanges::findOverlaps()]
#'
#' @return Either a SummarizedExperiment or GRanges object with additional
#' annotations provided by the supplied GRanges object.
#'
#' @examples
#' library(SummarizedExperiment)
#' data(rse_adar_ifn)
#' gr <- GRanges(rep(c("SSR3", "SPCS3"), c(5, 15)),
#'     IRanges(seq(1, 500, by = 25), width = 50),
#'     strand = "+"
#' )
#'
#' gr$feature <- sample(1:100, size = 20)
#' gr$id <- sample(LETTERS, size = 20)
#'
#' rse <- annot_from_gr(rse_adar_ifn, gr, c(feature_set = "feature", "id"))
#' rowData(rse)
#'
#' @importFrom S4Vectors aggregate unstrsplit
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<- seqlevels
#'
#' @export
annot_from_gr <- function(obj, gr, cols_to_map, RLE = TRUE, ...) {
    if (is(obj, "RangedSummarizedExperiment")) {
        gr_sites <- rowRanges(obj)
        return_se <- TRUE
    } else {
        gr_sites <- obj
        return_se <- FALSE
    }

    overlaps <- findOverlaps(gr_sites, gr, ...)

    if (!is.null(names(cols_to_map))) {
        names(cols_to_map) <- ifelse(names(cols_to_map) == "",
            cols_to_map,
            names(cols_to_map)
        )
    } else {
        names(cols_to_map) <- cols_to_map
    }

    for (i in seq_along(cols_to_map)) {
        col <- cols_to_map[[i]]
        col_id <- names(cols_to_map)[i]
        if (!col %in% names(mcols(gr))) {
            cli::cli_abort("{col} not present in mcols() of input")
        }
        mcols(gr)[[col]] <- as.character(mcols(gr)[[col]])
        x <- aggregate(gr,
            overlaps,
            tmp = unstrsplit(
                unique(eval(parse(text = col))),
                ","
            ),
            drop = FALSE
        )
        x$tmp <- ifelse(x$tmp == "", NA, x$tmp)
        mcols(gr_sites)[[col_id]] <- x$tmp
    }

    if(RLE) {
        col_nms <- names(cols_to_map)
        rls_clms <- lapply(col_nms, function(x) S4Vectors::Rle(mcols(gr_sites)[[x]]))
        mcols(gr_sites)[col_nms] <- rls_clms
    }

    if (return_se) {
        mcols(rowRanges(obj)) <- mcols(gr_sites)
    } else {
        mcols(obj) <- mcols(gr_sites)
    }

    obj
}
