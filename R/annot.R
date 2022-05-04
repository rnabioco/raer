#' @rdname annot_snps
#' @importFrom BSgenome snpsByOverlaps
#' @export
annot_snps.GRanges <- function(obj,
                       dbsnp,
                       chrom = NULL,
                       col_to_aggr = "RefSNP_id",
                       drop = FALSE){

  if(!is(dbsnp, "ODLT_SNPlocs")){
    stop(dbsnp, " not valid SNP package, please install SNP db")
  }

  if(!is.null(chrom)){
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
  mcols(sites)[col_to_aggr] <- aggregate(snps,
                            snp_overlaps,
                            snp = unstrsplit(eval(parse(text=col_to_aggr)),
                                             ","),
                            drop = FALSE)$snp

  seqlevelsStyle(sites) <- in_style

  if(drop){
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
                               drop = FALSE){
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
#'
#' @examples
#' example(create_se, echo = FALSE)
#' library(SummarizedExperiment)
#' gr <- GRanges(rep(c("SSR3", "SPCS3"), c(5, 15)),
#'              IRanges(seq(1, 500, by = 25), width=50),
#'              strand="+")
#' gr$feature <- sample(1:100, size = 20)
#' gr$id <- sample(LETTERS, size = 20)
#'
#' se <- annot_from_gr(se, gr, c(feature_set = "feature", "id"))
#' rowData(se)
#'
#' @importFrom S4Vectors aggregate unstrsplit
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @export
annot_from_gr <- function(obj, gr, cols_to_map){

  if(is(obj, "RangedSummarizedExperiment")){
    gr_sites <- rowRanges(obj)
    return_se <- TRUE
  } else {
    gr_sites <- obj
    return_se <- FALSE
  }

  missing_chroms <- setdiff(seqlevels(gr), seqlevels(gr_sites))
  if(length(missing_chroms) != 0){
    warning("The following chromosomes in gr are not present in obj\n",
            missing_chroms)
  }

  overlaps <- findOverlaps(gr_sites, gr, ignore.strand = TRUE)

  if(!is.null(names(cols_to_map))){
    names(cols_to_map) <- ifelse(names(cols_to_map) == "",
                                 cols_to_map,
                                 names(cols_to_map))
  } else {
    names(cols_to_map) <- cols_to_map
  }

  for(i in seq_along(cols_to_map)){
    col <- cols_to_map[[i]]
    col_id <- names(cols_to_map)[i]
    mcols(gr)[[col]] <- as.character(mcols(gr)[[col]])
    x <- aggregate(gr,
                   overlaps,
                   tmp = unstrsplit(eval(parse(text=col)),
                                    ","),
                   drop = FALSE)
    x$tmp <- ifelse(x$tmp == "", NA, x$tmp)
    mcols(gr_sites)[[col_id]] <- x$tmp
  }
  if(return_se){
    mcols(rowRanges(obj)) <- mcols(gr_sites)
  } else {
    mcols(obj) <- mcols(gr_sites)
  }

  obj
}
