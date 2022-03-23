
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
  mcols(sites) <- aggregate(snps,
                            snp_overlaps,
                            snp = unstrsplit(eval(parse(text=col_to_aggr)),
                                             ","),
                            drop = FALSE)

  seqlevelsStyle(sites) <- in_style

  if(drop){
    sites <- sites[!is.na(sites$snp)]
  }

  sites
}


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
#' to rowRanges of SummarizedExperiment. If multiple features overlap
#' then they will be concatenated as comma separated values.
#'
#' @param obj RangedSummarizedExperiment
#' @param gr GRanges with annotations to map to obj
#' @param cols_to_map list of columns from gr to map
#' to obj. If the list has names, the names will be used for
#' colnames in the obj
#'
#' @importFrom S4Vectors aggregate unstrsplit
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @export
annot_from_gr <- function(obj, gr, cols_to_map){

  gr_sites <- rowRanges(obj)
  site_style <- seqlevelsStyle(gr_sites)
  annot_style <- seqlevelsStyle(gr)

  if(site_style != annot_style){
    warning("gr and obj do not have the same chromosome convention\n",
            "attempting to coerce gr seqnames to obj")
    tryCatch( seqlevelsStyle(annot_style) <- annot_style,
              error = function(e) e,
              finally = print("unable to map chroms in gr to obj"))
  }

  overlaps <- findOverlaps(gr_sites, gr, ignore.strand = TRUE)

  if(is.null(names(cols_to_map))){
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
  mcols(rowRanges(obj)) <- mcols(gr_sites)
  obj
}
