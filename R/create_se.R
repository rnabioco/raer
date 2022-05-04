#' Create RangedSummarizedExperiment
#'
#' This function will take either a single result from running pileup_res or a list of
#' results (ie for different samples) from running pileup_res and will return a summarized
#' experiment object that contains assays for each column in the pileup_res output.
#' Currently this doesn't add any sample info and can only take a list of results.
#' Potential added functionality would be to allow the user to give an existing
#' SE object as input or to give sample meta data.
#'
#'
#' @param plps results from running get_pileup.R, can be one result, a list
#' of results, or a named list of results. If a named list is give, the colData
#' will be named using the names in the list.
#' @param rowdata_cols character vector of columns to store in rowData
#' . Values must be the same across all assays (excluding NA).
#' @param assay_cols character vector of columns to store as assays
#' @param sample_names OPTIONAL A list of names to be added to the SE object.
#' If no sample names are given and pileup_res is not a named list, then
#' default names (ie sample_1, sample_2, ..., sample_n) will be given and
#' a warning will be printed.
#'
#' @examples
#' library(SummarizedExperiment)
#' bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' bam2fn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#'
#' plps <- get_pileup(c(bamfn, bam2fn), fafn)
#' names(plps) <- c("sample1", "sample2")
#' se <- create_se(plps)
#'
#' assays(se)
#'
#' colData(se)
#'
#' rowRanges(se)
#'
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom IRanges extractList
#' @export

create_se <- function(plps,
                      rowdata_cols = c("Ref"),
                      assay_cols = c("Var", "nRef", "nVar", "nA", "nT", "nC", "nG"),
                      sample_names = NULL){
  if(!is.list(plps)){
    plps <- list(plps)
  }
  # Checks for sample names
  if(is.null(sample_names)){
    if(is.null(names(plps))){
      sample_names <- paste0("sample_", 1:length(plps))
    } else {
      sample_names <- names(plps)
    }
  } else {
    if(length(plps) != length(sample_names)){
        stop(paste0("You must provide the same number of sample names as pileup results!!!
                    You supplied ",
                    length(plps), " pileup results but supplied ",
                    length(sample_names), " sample names!!!"))

    }
  }

  # Find all ranges in the list of results
  all_ranges <- lapply(plps, function(x){
    mcols(x) <- NULL
    x
  })

  all_ranges <- GRangesList(all_ranges)
  all_ranges <- unique(unlist(all_ranges, use.names = FALSE))

  # Loop through all samples
  names(plps) <- sample_names
  se_list <- lapply(sample_names, function(sample_name){
    pileup_sample <- plps[[sample_name]]

    # Add NA where there is no data
    ranges_sample <- pileup_sample
    mcols(ranges_sample) <- NULL

    # Make a granges object of only the ranges that aren't in the current
    # sample
    hits <- findOverlaps(all_ranges, ranges_sample, type = "equal")
    grl <- extractList(ranges_sample, as(hits, "List"))
    unique_hits <- unlist(psetdiff(all_ranges, grl))

    # Add an empty metadata df to the newly made granges object
    meta_data <- matrix( NA,
                         nrow = NROW(unique_hits),
                         ncol = ncol(mcols(ranges_sample)))

    colnames(meta_data) <- colnames(mcols(ranges_sample))
    mcols(unique_hits) <- meta_data

    # Add the "empty" object to the
    pileup_sample <- c(pileup_sample, unique_hits)

    pileup_sample <- pileup_sample[order(match(pileup_sample, all_ranges))]

    assay_names <- colnames(mcols(pileup_sample))

    # Pull out data frame for each column in the results, these will become
    # the assays
    assay_list <- lapply(assay_names, function(x){
      return_df <- matrix(mcols(pileup_sample)[[x]])
      colnames(return_df) <- sample_name
      return_df
    })

    names(assay_list) <- assay_names

    rowRanges <- pileup_sample
    mcols(rowRanges) <- NULL
    colData <- data.frame(sample = sample_name)

    se <- SummarizedExperiment(assays = assay_list,
                               rowRanges = rowRanges, colData = colData)
    se
  })

  combined_se <- do.call(cbind, se_list)

  site_id <- paste0(seqnames(rowRanges(combined_se)),
                    "_",
                    start(rowRanges(combined_se)),
                    "_",
                    as.integer(as.factor(strand(rowRanges(combined_se)))))
  rowData(combined_se)$site_id <- site_id
  rownames(combined_se) <- site_id

  for(i in seq_along(rowdata_cols)){
    rc <- rowdata_cols[i]
    rowData(combined_se)[rc] <- apply(assay(combined_se, rc), 1, function(x) unique(x[!is.na(x)]))
    assay(combined_se, rc) <- NULL
  }
  assays_to_drop <- setdiff(assay_cols, names(assays(combined_se)))
  for(i in seq_along(assays_to_drop)) assay(combined_se,assays_to_drop[i]) <- NULL

  combined_se
}
