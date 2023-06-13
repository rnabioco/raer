#' @rdname annot_snps
#' @importFrom BSgenome snpsByOverlaps
#' @export
annot_snps.GRanges <- function(obj,
    dbsnp,
    chrom = NULL,
    col_to_aggr = "RefSNP_id",
    drop = FALSE,
    genome = NULL,
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
        snps$alt_alleles <- unstrsplit(snps$alt_alleles)

        snp_info <- aggregate(
            snps,
            snp_overlaps,
            snp = unstrsplit(eval(parse(text = col_to_aggr)), ","),
            ref_allele = unstrsplit(ref_allele, ","),
            alt_alleles = unstrsplit(alt_alleles, ","),
            drop = FALSE
        )
        snp_info$grouping <- NULL
        colnames(snp_info) <- c(col_to_aggr, "snp_ref_allele", "snp_alt_alleles")
        mcols(sites) <- cbind(mcols(sites), snp_info)
    }

    seqlevelsStyle(sites) <- in_style

    if (drop) {
        sites <- sites[!is.na(sites$snp)]
    }

    sites
}


#' @rdname annot_snps
#' @export
annot_snps.SummarizedExperiment <- function(obj,
    dbsnp,
    chrom = NULL,
    col_to_aggr = "RefSNP_id",
    drop = FALSE,
    genome = NULL,
    ...) {
    gr <- rowRanges(obj)
    res <- annot_snps.GRanges(gr,
        dbsnp,
        chrom = chrom,
        col_to_aggr = col_to_aggr,
        drop = drop,
        genome = genome
    )

    mcols(rowRanges(obj)) <- mcols(res)
    obj
}

#' Annotate a RangedSummarizedExperiment using Granges objects
#'
#' @description Utility function to map annotations from GRanges to rowData of
#'   SummarizedExperiment or GRanges object. If multiple features overlap then
#'   they will be concatenated as comma separated values.
#'
#' @param obj RangedSummarizedExperiment or GRanges object
#' @param gr GRanges with annotations to map to obj
#' @param cols_to_map character vector of columns from gr to map to obj. If the
#'   vector has names, the names will be the column names in the output obj
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
annot_from_gr <- function(obj, gr, cols_to_map, ...) {
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
    if (return_se) {
        mcols(rowRanges(obj)) <- mcols(gr_sites)
    } else {
        mcols(obj) <- mcols(gr_sites)
    }

    obj
}
