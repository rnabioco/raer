#' Create RangedSummarizedExperiment
#'
#' This function will take either a single result from running get_pileup() or a list of
#' results (ie for different samples) from running pileup_res and will return a summarized
#' experiment object that contains assays for specified columns in the get_pileup() output.
#'
#'
#' @param plps results from running [get_pileup()], can be one result, a list
#' of results, or a named list of results. If a named list is given, the colData
#' will be named using the names in the list.
#' @param assay_cols character vector of columns to store as assays
#' @param sample_names A list of names to be added to the SE object.
#' If no sample names are given and plps is not a named list, then
#' default names (ie sample_1, sample_2, ..., sample_n) will be given and
#' a warning will be printed.
#' @param sparse If TRUE, numeric matrices will be stored as sparseMatrices. All
#' missing values will be coerced to 0 in the sparseMatrices.
#' @param fill_na Numeric value to replace NAs in numeric matrices Should only be
#'  used when plps were computed independently with a min_nucleotide_count = 1,
#'  otherwise sites may be set to 0, although they may have coverage > 0 but less
#'  than the min_nucleotide_count parameter. Not applied when sparse = TRUE, which
#'  coerces missing values to 0.
#' @param verbose print information on progress
#'
#' @return
#' `RangedSummarizedExperiment` object populated with assays for each of the listed
#' `assay_cols`.
#' @examples
#' library(SummarizedExperiment)
#' bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' bam2fn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
#' fafn <- system.file("extdata", "human.fasta", package = "raer")
#'
#' bams <- rep(c(bamfn, bam2fn), each = 3)
#' sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
#' fp <- FilterParam(only_keep_variants = TRUE)
#' plps <- get_pileup(bams, fafn, filterParam = fp)
#' names(plps) <- sample_ids
#' se <- create_se(plps)
#' se$condition <- substr(se$sample, 1, 2)
#' assays(se)
#'
#' colData(se)
#'
#' rowRanges(se)
#'
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom IRanges extractList
#' @importFrom Matrix sparseMatrix
#' @export
create_se <- function(plps,
                      assay_cols = c("Var", "nRef", "nVar", "nA", "nT", "nC", "nG"),
                      sample_names = NULL,
                      sparse = FALSE,
                      fill_na = NULL,
                      verbose = FALSE) {
  if (!is.list(plps)) {
    plps <- list(plps)
  }
  # Checks for sample names
  if (is.null(sample_names)) {
    if (is.null(names(plps))) {
      sample_names <- paste0("sample_", 1:length(plps))
    } else {
      sample_names <- names(plps)
    }
  } else {
    if (length(plps) != length(sample_names)) {
      stop("You must provide the same number of sample names as pileup results",
        "You supplied ", length(plps), " pileup results but supplied ",
        length(sample_names), " sample names.")
    }
    names(plps) <- sample_names
  }

  stopifnot(all(assay_cols %in% colnames(mcols(plps[[1]]))))

  # determine all ranges for across plps,
  # store Ref nt in mcols (should be invariant)
  all_ranges <- unlist(as(plps, "GRangesList"), use.names = FALSE)
  ref_nt <- mcols(all_ranges)[["Ref"]]
  mcols(all_ranges) <- NULL
  mcols(all_ranges)$Ref <- ref_nt
  all_ranges <- sort(unique(all_ranges), ignore.strand = TRUE)
  assay_cols <- setdiff(assay_cols, "Ref")

  if (!sparse) {
    non_sparse_assays <- assay_cols
  } else {
    # character matrices aren't handled by Matrix
    ns <- !unlist(lapply(mcols(plps[[1]])[, assay_cols],
      function(x) is.numeric(x)))
    non_sparse_assays <- assay_cols[ns]
    if (any(ns)) {
      warning("sparseMatrices can only be generated for numeric values.\n",
        non_sparse_assays, " will be stored as dense matrices ",
        "which may consume excessive memory.")
    }
  }

  sparse_assays <- setdiff(assay_cols, non_sparse_assays)

  plp_assays <- fill_matrices(plps, non_sparse_assays, all_ranges, verbose)
  sp_plp_assays <- fill_sparse_matrices(plps, sparse_assays, all_ranges,
                                        use_hashmap = TRUE, verbose)
  plp_assays <- c(plp_assays, sp_plp_assays)

  se <- SummarizedExperiment(assays = plp_assays,
    rowRanges = all_ranges,
    colData = DataFrame(sample = sample_names))
  colnames(se) <- se$sample
  rownames(se) <- site_names(rowRanges(se))

  if (!is.null(fill_na) && !sparse) {
    assays(se) <- lapply(assays(se), function(x) {
      if (is.numeric(x)) x[is.na(x)] <- fill_na
      x
    })
  }
  se
}

# standardize and bind matrices from pileups
# into dense matrices
fill_matrices <- function(plps, assays, gr, verbose = TRUE) {
  if (length(assays) == 0) {
    return(NULL)
  }

  plps <- as(plps, "GRangesList")
  plp_assays <- lapply(assays, function(x) {
    matrix(NA,
      nrow = NROW(gr),
      ncol = length(plps))
  })
  names(plp_assays) <- assays

  for (i in seq_along(plps)) {
    if (verbose) {
      if (i %% 128 == 0) {
        message("processed ", i, " pileups")
      }
    }
    hits <- findOverlaps(gr, plps[[i]], type = "equal")
    hits <- queryHits(hits)
    for (j in seq_along(plp_assays)) {
      id <- names(plp_assays)[j]
      plp_assays[[j]][hits, i] <- mcols(plps[[i]])[[id]]
    }
  }
  plp_assays
}


bind_se <- function(ses){
  all_assays <- lapply(ses, function(x) names(assays(x)))
  assay_names <- unique(unlist(all_assays))
  stopifnot(all(unlist(lapply(all_assays,
                              function(x) {
                                length(setdiff(assay_names, x)) == 0
                                }))))

  assay_dat <- vector("list", length = length(assay_names))
  for(i in seq_along(assay_names)){
    an <- assay_names[i]
    an_assays <- lapply(ses, function(x) assay(x, an))
    assay_dat[[i]] <- cbind_sparse(an_assays)
    names(assay_dat)[i] <- an
  }

  cdat <- do.call(rbind, lapply(ses, colData))

  stopifnot(Reduce(identical, lapply(assay_dat, rownames)))
  rn <- rownames(assay_dat[[1]])
  rowrng <- unlist(GRangesList(lapply(ses, rowRanges)), use.names = FALSE)
  rowrng <- rowrng[rn, , drop = FALSE]
  SummarizedExperiment(assays = assay_dat,
                       rowRanges = rowrng,
                       colData = cdat)
}


# standardize and bind matrices from pileups
# into sparseMatrices
#' @import fastmap
fill_sparse_matrices <- function(plps, assays, gr, use_hashmap = TRUE, verbose = TRUE) {
  if (length(assays) == 0) {
    return(NULL)
  }
  plps <- as(plps, "GRangesList")

  plp_assays <- vector("list", length(assays))
  names(plp_assays) <- assays

  all_hits <- vector("list", length(plps))
  if (!use_hashmap) {
    ## use findOverlaps to find shared indexes
    for (i in seq_along(plps)) {
      if (verbose) {
        if (i %% 128 == 0) {
          message("processed ", i, " pileups")
        }
      }
      hits <- findOverlaps(gr, plps[[i]], type = "equal")
      all_hits[[i]] <- queryHits(hits)
    }
  } else {
    ## use rownames and hashmap to find shared indexes
    names(gr) <- site_names(gr)
    hm <- fastmap::fastmap()
    vals <- 1:length(gr)
    names(vals) <- names(gr)
    hm$mset(.list = as.list(vals))

    for (i in seq_along(plps)) {
      if (verbose) {
        if (i %% 128 == 0) {
          message("processed ", i, " pileups")
        }
      }
      if (length(plps[[i]]) == 0) {
        all_hits[[i]] <- integer(0)
      } else {
        all_hits[[i]] <- sort(unlist(hm$mget(site_names(plps[[i]])),
          use.names = FALSE))
      }
    }
    hm$reset()
    rm(hm)
  }

  for (i in seq_along(plp_assays)) {
    id <- names(plp_assays)[i]

    # drop zero value entries, otherwise stored in sparseMatrix
    vals <- lapply(plps, function(x) mcols(x)[, id])
    to_keep <- lapply(vals, function(v) v != 0)
    hits <- all_hits
    for (vi in seq_along(vals)) {
      vals[[vi]] <- vals[[vi]][to_keep[[vi]]]
      hits[[vi]] <- hits[[vi]][to_keep[[vi]]]
    }
    x <- cpp_fill_sparse_matrix(vals, hits)
    plp_assays[[i]] <- sparseMatrix(i = x[, 1], j = x[, 2], x = x[, 3],
      dims = c(length(gr), length(plps)))
  }
  plp_assays
}


# adapted from user6376297
# https://stackoverflow.com/a/56547070/6276041
# cbind sparse matrices with differing row entries
cbind_sparse <- function(mats){
  stopifnot(all(unlist(lapply(mats, function(x) is(x, "dgCMatrix")))))
  rn <- unique(unlist(lapply(mats, rownames)))
  cn <- unique(unlist(lapply(mats, colnames)))
  nvals <- sum(unlist(lapply(mats, Matrix::nnzero)))
  x <- integer(nvals)
  i <- integer(nvals)
  j <- integer(nvals)

  init_x <- 1
  for (idx in seq_along(mats)) {
    cindnew <- match(colnames(mats[[idx]]),cn)
    rindnew <- match(rownames(mats[[idx]]),rn)
    ind <- summary(mats[[idx]])
    nval <- nrow(ind)
    end <- nval + init_x - 1
    i[init_x:end] <- rindnew[ind$i]
    j[init_x:end] <- cindnew[ind$j]
    x[init_x:end] <- ind$x
    init_x <- end + 1
  }

  sparseMatrix(i = i,
               j = j,
               x = x,
               dims=c(length(rn), length(cn)),
               dimnames=list(rn, cn))
}


#' @importFrom stringr str_c
site_names <- function(gr) {
  if (length(gr) == 0) {
    return(NULL)
  }
  stringr::str_c(seqnames(gr),
    "_",
    start(gr),
    "_",
    as.integer(strand(gr)))
}
