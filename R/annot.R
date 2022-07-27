#' @rdname annot_snps
#' @importFrom BSgenome snpsByOverlaps
#' @export
annot_snps.GRanges <- function(obj,
                               dbsnp,
                               chrom = NULL,
                               col_to_aggr = "RefSNP_id",
                               drop = FALSE,
                               ...) {

  if (!is(dbsnp, "ODLT_SNPlocs")) {
    stop(dbsnp, " not valid SNP package, please install SNP db")
  }

  if (!is.null(chrom)) {
    sites <- obj[seqnames(obj) == chrom]
  } else {
    sites <- obj
  }

  in_style <- seqlevelsStyle(sites)
  snp_style <- seqlevelsStyle(dbsnp)
  seqlevelsStyle(sites) <- snp_style

  # returns GPos for each snp
  snps <- snpsByOverlaps(dbsnp, sites)

  # now annot if site is a SNP
  snp_overlaps <- findOverlaps(sites, snps, ignore.strand = TRUE)
  if (length(snp_overlaps) == 0) {
    mcols(sites)[col_to_aggr] <- NA
    seqlevelsStyle(sites) <- in_style
    return(sites)
  }
  mcols(sites)[col_to_aggr] <- aggregate(snps,
    snp_overlaps,
    snp = unstrsplit(eval(parse(text = col_to_aggr)),
      ","),
    drop = FALSE)$snp

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
                                            ...) {
  gr <- rowRanges(obj)
  res <- annot_snps.GRanges(gr,
    dbsnp,
    chrom = chrom,
    col_to_aggr = col_to_aggr,
    drop = drop)

  res <- cbind(mcols(gr), mcols(res))
  mcols(rowRanges(obj)) <- res
  obj
}

#' Annotate a RangedSummarizedExperiment using Granges objects
#'
#' @description Utility function to map annotations from GRanges
#' to rowData of SummarizedExperiment or GRanges object . If multiple features
#' overlap then they will be concatenated as comma separated values.
#'
#' @param obj RangedSummarizedExperiment or GRanges object
#' @param gr GRanges with annotations to map to obj
#' @param cols_to_map character vector of columns from gr to map
#' to obj. If the vector has names, the names will be the
#' column names in the output obj
#' @param ... additional arguments to pass to [GenomicRanges::findOverlaps()]
#'
#' @examples
#' example(create_se, echo = FALSE)
#' library(SummarizedExperiment)
#' gr <- GRanges(rep(c("SSR3", "SPCS3"), c(5, 15)),
#'   IRanges(seq(1, 500, by = 25), width = 50),
#'   strand = "+")
#' gr$feature <- sample(1:100, size = 20)
#' gr$id <- sample(LETTERS, size = 20)
#'
#' se <- annot_from_gr(se, gr, c(feature_set = "feature", "id"))
#' rowData(se)
#'
#' @importFrom S4Vectors aggregate unstrsplit
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<- seqlevels
#' @export
annot_from_gr <- function(obj, gr, cols_to_map, ...) {

  if (is(obj, "RangedSummarizedExperiment")) {
    gr_sites <- rowRanges(obj)
    return_se <- TRUE
  } else {
    gr_sites <- obj
    return_se <- FALSE
  }

  # missing_chroms <- setdiff(seqlevels(gr), seqlevels(gr_sites))
  # if(length(missing_chroms) != 0){
  #   warning("The following chromosomes in gr are not present in obj\n",
  #           missing_chroms)
  # }

  overlaps <- findOverlaps(gr_sites, gr, ...)

  if (!is.null(names(cols_to_map))) {
    names(cols_to_map) <- ifelse(names(cols_to_map) == "",
      cols_to_map,
      names(cols_to_map))
  } else {
    names(cols_to_map) <- cols_to_map
  }

  for (i in seq_along(cols_to_map)) {
    col <- cols_to_map[[i]]
    col_id <- names(cols_to_map)[i]
    if (!col %in% names(mcols(gr))) {
      stop(col, " not present in mcols() of input")
    }
    mcols(gr)[[col]] <- as.character(mcols(gr)[[col]])
    x <- aggregate(gr,
      overlaps,
      tmp = unstrsplit(eval(parse(text = col)),
        ","),
      drop = FALSE)
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

#' @importFrom BiocGenerics unlist
annot_sites <- function(obj,
                         txdb,
                         features = c("gene",
                                      "transcript",
                                      "three_utr",
                                      "five_utr",
                                      "exon",
                                      "cds",
                                      "intron")) {
  if (is(obj, "RangedSummarizedExperiment")) {
    gr_sites <- rowRanges(obj)
    return_se <- TRUE
  } else {
    gr_sites <- obj
    return_se <- FALSE
  }

  txfeatures <- get_txdb_features(txdb, features)
  annots <- lapply(seq_along(txfeatures), function(i){
    x <- txfeatures[[i]]
    id <- names(txfeatures)[i]
    if(id == "gene"){
      id_col <- "gene_id"
      gr <- suppressWarnings(annot_from_gr(gr_sites, x, id_col))
      res <- DataFrame(gene_id = gr$gene_id,
                       region = ifelse(is.na(mcols(gr)[[id_col]]), NA, id))
    } else {
      id_col <- "tx_name"
      gr <- suppressWarnings(annot_from_gr(gr_sites, x, id_col))
      res <- DataFrame(region = ifelse(is.na(mcols(gr)[[id_col]]), NA, id))
    }
    res
  })
  names(annots) <- features

  annot_regions <- do.call(cbind, annots)
  regions <- paste_cols(annot_regions)

  if("gene" %in% features){
    annot_cols <- annots$gene
    annot_cols$region <- regions
  } else {
    annot_cols <- DataFrame(region = regions)
  }
  mcols(gr_sites) <- cbind(mcols(gr_sites), annot_cols)

}

# templated on approach from RCAS::getTxdbFeaturesFromGRanges
get_txdb_features <- function(txdb, features){
  txfeatures <- list(
    gene = NULL,
    transcript = NULL,
    exon = NULL,
    intron = NULL,
    cds = NULL,
    three_utr = NULL,
    five_utr = NULL)

  features_requested <- setdiff(features, names(txfeatures))
  if(length(features_requested) > 0){
    input_msg <- paste(features_requested, collapse = ", ")
    req_msg <- paste(names(txfeatures), collapse = ", ")
    stop("features ", input_msg," must be one of ", req_msg)
  }
  if("transcript" %in% features){
    txfeatures$transcript <- GenomicFeatures::transcripts(txdb)
  }

  if("gene" %in% features) {
    tmp <- GenomicFeatures::genes(txdb, single.strand.genes.only=FALSE)
    gene <- BiocGenerics::unlist(tmp)
    gene$gene_id <- names(gene)
    txfeatures$gene <- unname(gene)
  }

  if("exon" %in% features) {
    tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
    exons <- BiocGenerics::unlist(tmp)
    exons$tx_name <- names(exons)
    txfeatures$exon <- unname(exons)
  }

  if("intron" %in% features) {
    tmp <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
    introns <- BiocGenerics::unlist(tmp)
    introns$tx_name <- names(introns)
    txfeatures$intron <- unname(introns)
  }

  if("five_utr" %in% features) {
    tmp <- range(GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
    fiveUTRs <- BiocGenerics::unlist(tmp)
    fiveUTRs$tx_name <- names(fiveUTRs)
    txfeatures$five_utr <- unname(fiveUTRs)
  }

  if("three_utr" %in% features){
    tmp <- range(GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE))
    threeUTRs <- BiocGenerics::unlist(tmp)
    threeUTRs$tx_name <- names(threeUTRs)
    txfeatures$three_utr <- unname(threeUTRs)
  }

  if("cds" %in% features){
    tmp <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
    cds <- BiocGenerics::unlist(tmp)
    cds$tx_name <- names(cds)
    txfeatures$cds <- unname(cds)
  }

  txfeatures[!unlist(lapply(txfeatures, is.null))]
}


# paste strings in multiple columns into a chracter vector

paste_cols <- function(dframe, cols = NULL, sep = ",", na.rm = TRUE){
  if(is.null(cols)) {
    cols <- colnames(dframe)
  }
  string_mat <- t(unname(sapply(dframe[cols], as.character)))
  if(na.rm){
    res <- apply(string_mat, 2, function(x) paste(x[!is.na(x)], collapse = sep))
  } else {
    res <- apply(string_mat, 2, function(x) paste(x, collapse = sep))
  }
  res
}
