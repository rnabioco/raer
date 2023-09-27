#' Adds editing frequencies
#'
#' @description Adds editing frequencies to an existing
#' [RangedSummarizedExperiment] object (created by [pileup_sites()]). The
#' [RangedSummarizedExperiment] with a new assay for editing frequencies
#' for each site (`edit_freq`), depth of coverage computed
#' using the indicated edited nucleotides (`depth`) and new `colData`
#' columns with the number of edited sites (`n_sites`) and the
#' fraction of edits (`edit_idx`) is returned.
#'
#' @param rse A [RangedSummarizedExperiment] object created by [pileup_sites()]
#' @param edit_from This should correspond to a nucleotide or assay
#'  (`A`, `C`, `G`, `T`, `Ref`, or `Alt`) you expect in the reference. Ex. for A to I
#'   editing events, this would be `A`.
#' @param edit_to This should correspond to a nucleotide or assay
#'  (`A`, `C`, `G`, `T`, `Ref`, or `Alt`)  you expect in the editing site. Ex. for A
#'   to I editing events, this would be `G`.
#' @param drop If `TRUE`, the [RangedSummarizedExperiment] returned will only retain sites
#'   matching the specified `edit_from` and `edit_to` bases.
#' @param replace_na If `TRUE`, `NA` and `NaN` editing frequencies will be coerced to
#'   `0`.
#' @param edit_frequency  The edit frequency cutoff used when calculating the
#'   number of sites. Set to `0` to require any non-zero editing frequency. The
#'   number of sites is stored as `n_sites` in the `colData`.
#' @param min_count The minimum number of reads required when enumerating number
#'   of editing sites detected.
#'
#' @return
#' [RangedSummarizedExperiment] supplemented with `edit_freq` and `depth` assay.
#'
#' @examples
#' library(SummarizedExperiment)
#' rse_adar_ifn <- mock_rse()
#' rse <- calc_edit_frequency(rse_adar_ifn)
#' assay(rse, "edit_freq")[1:5, ]
#'
#' @import SummarizedExperiment
#' @importFrom Matrix colSums
#' @export
calc_edit_frequency <- function(rse,
    edit_from = "A",
    edit_to = "G",
    drop = FALSE,
    replace_na = TRUE,
    edit_frequency = 0,
    min_count = 1) {

    valid_assays <- c("A", "T", "C", "G", "Ref", "Alt")

    if (!(edit_from %in% valid_assays) | !(edit_to %in% valid_assays)) {
        cli::cli_abort("`edit_to` and `edit_from` must be one of: ",
                       "'A', 'T', 'C', 'G', 'Ref', or 'Alt'")
    }

    from_col <- paste0("n", edit_from)
    to_col <- paste0("n", edit_to)

    if (drop && from_col != "nRef") {
        rse <- rse[mcols(rowRanges(rse))$REF == edit_from, ]
    }

    if ("depth" %in% names(assays(rse))) {
        cli::cli_alert_info(
            "depth has been overwritten with sum of {to_col} and {from_col} assays"
            )
    }

    assay(rse, "depth") <- assay(rse, to_col) + assay(rse, from_col)
    no_depth <- Matrix::rowSums(assay(rse, "depth")) == 0
    if (any(no_depth)) {
        cli::cli_alert_info(
            "{sum(no_depth)} sites had no coverage for calculating editing\n",
            "    these sites have been removed"
        )
        rse <- rse[!no_depth, ]
    }

    if (is(assay(rse, to_col), "sparseMatrix") ||
        is(assay(rse, from_col), "sparseMatrix")) {
        if (replace_na) {
            cli::cli_abort("NA values cannot be stored in sparseMatrices")
        }
        # compute editing frequencies, only at non-zero depth positions
        # zero depth positions are not stored in matrix
        # if coerced to simple matrix these will have editing
        # frequencies of 0
        idx <- Matrix::which(assay(rse, "depth") > 0, arr.ind = TRUE)
        res <- Matrix::sparseMatrix(idx[, 1],
            idx[, 2],
            x = assay(rse, to_col)[idx] /
                assay(rse, "depth")[idx]
        )
        dimnames(res) <- dimnames(assay(rse, from_col))
    } else {
        res <- assay(rse, to_col) / assay(rse, "depth")
        if (replace_na) {
            res[is.na(res)] <- 0
        }
    }
    assay(rse, "edit_freq") <- res
    rse <- count_edits(rse, edit_frequency, min_count, edit_from, edit_to)
    rse
}

#' Counts edits
#'
#' @description Counts edits per sample and add new colData columns with the
#'   number of edited sites (n_sites) and the  fraction of edits (edit_idx).
#'   This function should be called by `calc_edit_frequency` and is not meant to
#'   be used directly.
#'
#' @param se A SummarizedExperiment object created by `merge_pileups` and
#'   processed by `calc_edit_frequency`
#' @param edit_from OPTIONAL if not using a pre-built type, you can specify your
#'   own editing. This should be a nucleotide (A, C, G, or T) and should
#'   correspond to the nucleotide you expect in the reference. Ex. for A to I
#'   editing events, this would be "A". If type is not "AI", both edit from and
#'   edit_to must be set.
#' @param edit_to OPTIONAL if not using a pre-built type, you can specify your
#'   own editing. This should be a nucleotide (A, C, G, or T) and should
#'   correspond to the nucleotide you expect after the editing event. Ex. for A
#'   to I editing events, this would be "G". If type is not "AI", both edit from
#'   and edit_to must be set.
#' @param edit_frequency OPTIONAL the edit frequency used to determine the
#'   number of sites. Default is 0.01.
#' @param min_count OPTIONAL the number of reads used to determine the number of
#'   edited sites. Default is 10.
#'
#' @import SummarizedExperiment
#' @importFrom Matrix colSums
#' @noRd
#' @keywords internal
count_edits <- function(se, edit_frequency = 0.01, min_count = 10,
    edit_from = NULL, edit_to = NULL) {
    n_pass_filter <- Matrix::colSums((assay(se, "edit_freq") > edit_frequency) &
        ((assay(se, paste0("n", edit_from)) +
            assay(se, paste0("n", edit_to))) >=
            min_count))

    colData(se)$n_sites <- n_pass_filter

    edit_idx <- Matrix::colSums(assay(se, paste0("n", edit_to))) /
        (Matrix::colSums(assay(se, paste0("n", edit_from))) +
            Matrix::colSums(assay(se, paste0("n", edit_to))))

    colData(se)$edit_idx <- edit_idx

    se
}

#' Make summarized experiment object for differential editing analysis
#'
#' @description Generates a [RangedSummarizedExperiment] object for use with
#' `edgeR` or `DESeq2` . Will generate a `counts` assay with
#' a matrix formatted with 2 columns per sample,
#' representing the reference and editing allele counts.
#'
#' @param rse A [RangedSummarizedExperiment] object
#' @param edit_from This should correspond to a nucleotide or assay
#'  (`A`, `C`, `G`, `T`, `Ref`, or `Alt`) you expect in the reference. Ex. for A to I
#'   editing events, this would be `A`.
#' @param edit_to This should correspond to a nucleotide or assay
#'  (`A`, `C`, `G`, `T`, `Ref`, or `Alt`) you expect in the editing site. Ex. for A
#'   to I editing events, this would be `G`.
#' @param min_prop The minimum required proportion of reads edited at a site. At
#'   least `min_samples` need to pass this to keep the site.
#' @param max_prop The maximum allowable proportion of reads edited at a site. At
#'   least `min_samples` need to pass this to keep the site.
#' @param min_samples The minimum number of samples passing the `min_prop` and
#'   `max_prop` cutoffs to keep a site.
#'
#' @import SummarizedExperiment
#'
#' @returns  [RangedSummarizedExperiment] for use with `edgeR` or
#'  `DESeq2`. Contains a `counts` assay with a matrix formatted
#'  with 2 columns per sample (ref and alt counts).
#'
#' @examples
#' library(SummarizedExperiment)
#' rse_adar_ifn <- mock_rse()
#' rse <- calc_edit_frequency(rse_adar_ifn)
#' dse <- make_de_object(rse, min_samples = 1)
#' assay(dse, "counts")[1:5, ]
#' dse
#' @export
make_de_object <- function(rse,
    edit_from = "A",
    edit_to = "G",
    min_prop = 0.0,
    max_prop = 1.0,
    min_samples = 1) {

    valid_assays <- c("A", "T", "C", "G", "Ref", "Alt")

    if (!(edit_from %in% valid_assays) | !(edit_to %in% valid_assays)) {
        cli::cli_abort("`edit_to` and `edit_from` must be one of: ",
                       "'A', 'T', 'C', 'G', 'Ref', or 'Alt'")
    }

    # Only keep locations that pass cutoffs in a certain number of samples
    if(!"edit_freq" %in% assayNames(rse)) {
        cli::cli_abort("`edit_freq` not present in assay, please run calc_edit_frequency()")
    }
    pass_cutoff <- (assay(rse, "edit_freq") >= min_prop) &
        (assay(rse, "edit_freq") <= max_prop)
    rse <- rse[rowSums(pass_cutoff) >= min_samples, ]

    # Set the ref and alternate allele and create a count table with both
    ref <- assay(rse, paste0("n", edit_from))
    colnames(ref) <- paste0(colnames(ref), "_ref")
    alt <- assay(rse, paste0("n", edit_to))
    colnames(alt) <- paste0(colnames(alt), "_alt")
    res <- cbind(ref, alt)
    mdata <- colData(rse)

    # Join the meta data for all samples
    ref_mdata <- alt_mdata <- mdata
    rownames(ref_mdata) <- colnames(ref)
    ref_mdata$count <- "ref"
    rownames(alt_mdata) <- colnames(alt)
    alt_mdata$count <- "alt"
    mdata <- rbind(ref_mdata, alt_mdata)
    mdata$count <- factor(mdata$count)

    # create a new SummarizedExperiment
    res <- SummarizedExperiment(
        assays = list(counts = res),
        colData = mdata
    )
    rowRanges(res) <- rowRanges(rse)
    res
}

#' Perform differential editing
#'
#' @description Use `edgeR` or `DESeq2` to perform differential editing
#'   analysis. This will work for simple designs that have 1 treatment and 1
#'   control. For more complex designs, we suggest you perform your own.
#'
#'   At the moment, this function will only find editing events specific to the
#'   treatment.
#'
#' @param deobj A [RangedSummarizedExperiment] object prepared for differential
#' editing analysis by [make_de_object()]
#' @param test Indicate if `edgeR` or `DESeq2` should be run.
#' @param sample_col The name of the column from `colData(deobj)` that
#'   contains your sample information. Default is sample. If you do not have a
#'   column named "sample", you must provide the appropriate sample column
#' @param condition_col The name of the column from `colData(deobj)` that
#'   contains your treatment information. Default is condition, If you do not
#'   have a column named "condition", you must provide the appropriate condition
#'   column
#' @param condition_control The name of the control condition. This must be a
#'   variable in your `condition_col` of `colData(deobj)`. No default provided.
#' @param condition_treatment The name of the treatment condition. This must be
#'   a variable in your `condition_col` of `colData(deobj)`.
#'
#' @examples
#' library(SummarizedExperiment)
#' bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
#' bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
#' fafn <- raer_example("human.fasta")
#'
#' bams <- rep(c(bamfn, bam2fn), each = 3)
#' sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
#' names(bams) <- sample_ids
#'
#' fp <- FilterParam(only_keep_variants = TRUE)
#' rse <- pileup_sites(bams, fafn, param = fp)
#' rse$condition <- substr(rse$sample, 1, 2)
#'
#' rse <- calc_edit_frequency(rse)
#' dse <- make_de_object(rse)
#' res <- find_de_sites(dse, condition_control = "WT", condition_treatment = "KO")
#' res$sig_results[1:3, ]
#'
#' @returns A named list:
#'   - `de_obj`: The `edgeR` or `deseq` object used for differential editing analysis
#'   - `results_full`: Unfiltered differential editing results
#'   - `sig_results`: Filtered differential editing (FDR < 0.05)
#'   - `model_matrix`: The model matrix used for generating DE results
#'
#' @importFrom stats model.matrix
#' @export
find_de_sites <- function(deobj,
                       test = c("edgeR", "DESeq2"),
                       sample_col = "sample",
                       condition_col = "condition",
                       condition_control = NULL,
                       condition_treatment = NULL) {
    # Make sure all variables are present
    if (!sample_col %in% colnames(colData(deobj))) {
        cli::cli_abort(
            c(
                "somple_col must be a column in the colDat of your deobj.",
                "'{sample_col}' not found in colnames(colData(deobj))"
            )
        )
    }
    if (!condition_col %in% colnames(colData(deobj))) {
        cli::cli_abort(
            c(
                "condition_col must be a column in the colData of your deobj.",
                "'{condition_col}' not found in colnames(colData(deobj))!"
            )
        )
    }

    # Rename columns based on the input
    if (sample_col != "sample") {
        colData(deobj)$sample <- NULL
    }
    if (condition_col != "condition") {
        colData(deobj)$condition <- NULL
    }

    new_columns <- as.data.frame(colData(deobj))
    names(new_columns)[names(new_columns) == condition_col] <-
        "condition"
    names(new_columns)[names(new_columns) == sample_col] <- "sample"
    new_columns <- new_columns[c("sample", "condition")]

    colData(deobj) <- cbind(colData(deobj), new_columns)

    # Check that condition_control and condition_treatment are correct
    if (is.null(condition_control)) {
        cond_options <-
            paste(unique(colData(deobj)$condition), collapse = ", ")
        cli::cli_abort(
            c(
                "condition_control must be set. This should be the level of",
                " your meta data that corresponds to your control. Possible",
                " options from your experiment are: {cond_options}"
            )
        )
    }
    if (is.null(condition_treatment)) {
        cond_options <-
            paste(unique(colData(deobj)$condition), collapse = ", ")
        cli::cli_abort(
            c(
                "condition_treatment must be set. This should be the level of",
                " your meta data that corresponds to your control. Possible",
                " options from your experiment are: {cond_options}"
            )
        )
    }

    # Check that the treatment and control are in the object
    if (!condition_control %in% colData(deobj)$condition) {
        cond_options <-
            paste(unique(colData(deobj)$condition), collapse = ", ")
        cli::cli_abort(
            c(
                "condition_control must be a column in your deobj colData.",
                "'{condition_control}' not found in the levels of the condition",
                " column of colData(deobj). Possible",
                " options from your experiment are: {cond_options}"
            )
        )
    }
    if (!condition_treatment %in% colData(deobj)$condition) {
        cond_options <-
            paste(unique(colData(deobj)$condition), collapse = ", ")
        cli::cli_abort(
            c(
                "condition_treatment must be a column in your deobj colData ",
                "'{condition_treatment}' not found in the levels of the condition",
                " column of colData(deobj)! Possible",
                " options from your experiment are: {cond_options}"
            )
        )
    }
    test <- match.arg(test)
    if (test == "edgeR") {
        de_fun <- run_edger
    } else if (test == "DESeq2") {
        de_fun <- run_deseq2
    } else {
        cli::cli_abort(c(
            "Unrecognized type: '{test}'. type must be either edgeR",
            " or DESeq2."
        ))
    }
    results <- de_fun(deobj, condition_control, condition_treatment)
    results
}

# Perform differential editing with DESeq2
#
# @description Uses DESeq2 to perform differential editing analysis. This will
#   work for simple designs that have 1 treatment and 1 control. For more
#   complex designs, we suggest you perform your own. It will test if your
#   sample column makes the model matrix not full rank. If that happens, the
#   model matrix will be modified to be full rank. This is not intended to be
#   called directly by the user, instead, this should be called by `find_de_sites`
#
#   At the moment, this function will only find editing events specific to the
#   treatment, but it will be pretty straight forward to add other possible
#   return values.
#
# @param deobj A SummarizedExperiment object prepared for de by `make_de_object`
# @param condition_control The name of the control condition. This must be a
#   variable in your condition_col of colData(deobj). No default provided.
# @param condition_treatment The name of the treatment condition. This must be
#   a variable in your condition_col of colData(deobj).
#
#
run_deseq2 <- function(deobj, condition_control = NULL,
    condition_treatment = NULL) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        cli::cli_abort("Package \"DESeq2\" needed for differential editing.")
    }

    design <- ~ 0 + condition:sample + condition:count

    # See if the design is full rank, if not, remove sample info
    test_mat <- try(DESeq2::DESeqDataSetFromMatrix(
        countData = assay(
            deobj,
            "counts"
        ),
        colData = colData(deobj),
        design = design
    ), silent = TRUE)

    if (is(test_mat, "try-error")) {
        sample <- deobj$sample
        condition <- deobj$condition
        count <- deobj$count

        design <- model.matrix(~ 0 + sample + condition:count)
        design <- design[, !grepl("countref", colnames(design))]
        mod_mat <- design
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = assay(deobj, "counts"),
        colData = colData(deobj),
        design = design
    )

    # We don't want size factors because we are
    # looking at ratios within a sample
    DESeq2::sizeFactors(dds) <- rep(1, nrow(colData(deobj)))

    dds <- DESeq2::DESeq(dds)

    if (!is(test_mat, "try-error")) {
        mod_mat <- model.matrix(design(dds), colData(dds))
    }

    # Pull out the model matrix for all comparisons of interest
    alt_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
        dds$count == "alt", ])
    ref_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
        dds$count == "ref", ])

    alt_control <- colMeans(mod_mat[dds$condition == condition_control &
        dds$count == "alt", ])
    ref_control <- colMeans(mod_mat[dds$condition == condition_control &
        dds$count == "ref", ])

    # This finds editing specific to the treatment condition
    treatment_vs_control <- DESeq2::results(dds,
        contrast = (alt_treatment -
            ref_treatment) -
            (alt_control - ref_control)
    )

    deseq_res <- as.data.frame(treatment_vs_control)
    deseq_res <- deseq_res[deseq_res$padj < 0.05, ]
    deseq_res <- deseq_res[order(deseq_res$log2FoldChange, decreasing = TRUE), ]

    return(list(
        de_obj = dds,
        results_full = treatment_vs_control,
        sig_results = deseq_res,
        model_matrix = mod_mat
    ))
}

# Perform differential editing with edgeR
#
# @description Uses edgeR to perform differential editing analysis.
# This will work for simple designs that have 1 treatment and 1 control. For
# more complex designs, we suggest you perform your own. It will test if your
# sample column makes the
# model matrix not full rank. If that happens, the model matrix will be
# modified to be full rank. This is not intended to be called directly by the
# user, instead, this should be called by `calc_differential_editing`
#
# At the moment, this function will only find editing events specific to the
# treatment, but it will be pretty straight forward to add other possible
# return values.
#
# @param deobj A SummarizedExperiment object prepared for de by `make_de_object`
# @param condition_control The name of the control condition. This must be a
#   variable in your condition_col of colData(deobj). No default provided.
# @param condition_treatment The name of the treatment condition. This must be
#   a variable in your condition_col of colData(deobj).
run_edger <- function(deobj, condition_control = NULL,
    condition_treatment = NULL) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        cli::cli_abort("Package \"edgeR\" needed to run differential analysis.")
    }

    sample <- deobj$sample
    condition <- deobj$condition
    count <- deobj$count

    design <- model.matrix(~ 0 + condition:sample + condition:count)

    # Check if full rank, if not, fix
    if (!limma::is.fullrank(design)) {
        design <- model.matrix(~ 0 + sample + condition:count)
        design <- design[, !grepl("countref", colnames(design))]
    }

    dge <- edgeR::DGEList(assay(deobj, "counts"),
                          lib.size = rep(1, ncol(deobj)))
    dge <- edgeR::estimateDisp(dge)
    fit <- edgeR::glmFit(dge, design)

    # Pull out the model matrix for all comparisons of interest
    alt_treatment <- colMeans(design[deobj$condition == condition_treatment &
        deobj$count == "alt", ])
    ref_treatment <- colMeans(design[deobj$condition == condition_treatment &
        deobj$count == "ref", ])

    alt_control <- colMeans(design[deobj$condition == condition_control &
        deobj$count == "alt", ])
    ref_control <- colMeans(design[deobj$condition == condition_control &
        deobj$count == "ref", ])

    treatment_vs_control <- edgeR::glmLRT(fit,
                                          contrast = (
                                              alt_treatment - ref_treatment) -
                                              (alt_control - ref_control)
                                          )

    treatment_vs_control <- edgeR::topTags(treatment_vs_control,
        n = nrow(treatment_vs_control)
    )

    edger_res <- as.data.frame(treatment_vs_control)
    edger_res <- edger_res[edger_res$FDR < 0.05, ]
    edger_res <- edger_res[order(edger_res$PValue), ]

    return(list(
        de_obj = fit,
        results_full = treatment_vs_control,
        sig_results = edger_res,
        model_matrix = design
    ))
}

#' Identify sites with differential editing between cells in single cell datasets
#'
#' @description
#' Compare editing frequencies between clusters or celltypes. REF and ALT counts
#' from each cluster are pooled to create pseudobulk estimates. Each pair of clusters
#' are compared using fisher exact tests. Statistics are aggregated across each pairwise
#' comparison using [scran::combineMarkers].
#'
#' @param sce [SingleCellExperiment] object with `nRef` and `nAlt` assays.
#' @param group column name from colData used to define groups to compare.
#' @param BPPARAM BiocParallel backend for control how paralllel computations are
#' performed.
#' @param ... Additional arguments passed to [scran::combineMarkers]
#'
#' @returns
#' A named list of [DataFrame]s containing results for each cluster specified by `group`.
#' The difference in editing frequencies between cluster pairs are denoted as `dEF`.
#' See [scran::combineMarkers] for a description of additional output fields.
#'
#'
#'
#' @examples
#'
#' ### generate example data ###
#'
#' library(Rsamtools)
#' library(GenomicRanges)
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#'
#' gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
#' gr$REF <- c(rep("A", 4), "T")
#' gr$ALT <- c(rep("G", 4), "C")
#'
#' cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
#' cbs <- na.omit(cbs)
#'
#' outdir <- tempdir()
#' bai <- indexBam(bam_fn)
#'
#' fp <- FilterParam(library_type = "fr-second-strand")
#' sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
#'
#' # mock some clusters
#' set.seed(42)
#' sce$clusters <- paste0("cluster_", sample(1:3, ncol(sce), replace = TRUE))
#' res <- find_scde_sites(sce, "clusters")
#' res[[1]]
#' @export
find_scde_sites <- function(
        sce,
        group,
        BPPARAM = SerialParam(),
        ...) {

    if (!requireNamespace("scran", quietly = TRUE)) {
        cli::cli_abort("Package \"scran\" needed for differential editing.")
    }

    if (!requireNamespace("scuttle", quietly = TRUE)) {
        cli::cli_abort("Package \"scran\" needed for differential editing.")
    }

    if(!is(sce, "SingleCellExperiment")) {
        cli::cli_abort("sce must be a SingleCellExperiment object")
    }

    if(!all(c("nRef", "nAlt") %in% assayNames(sce))) {
        cli::cli_abort("sce must contain nRef and nAlt assays")
    }

    if(length(group) != 1) {
        cli::cli_abort("group must be a single value")
    }

    if(!(group %in% colnames(colData(sce)))){
        cli::cli_abort("{group} not found in colData")
    }

    if(any(is.na(sce[[group]]))) {
        cli::cli_abort("NA values not allowed in group column")
    }

    if(!"depth" %in% names(assays(sce))) {
        assay(sce, "depth") <- assay(sce, "nRef") + assay(sce, "nAlt")
    }

    no_depth_cells <-  colSums(assay(sce, "depth")) == 0
    if(sum(no_depth_cells) > 0) {
        cli::cli_alert(c("{sum(no_depth_cells)} cells had no REF or ALT counts ",
                         "and were excluded from the analysis"))
        sce <- sce[, !no_depth_cells]
    }

    assay(sce, "depth") <- scuttle::normalizeCounts(sce,
                                                    assay.type = "depth",
                                                    BPPARAM = BPPARAM)

    nref <- scuttle::summarizeAssayByGroup(sce,
                                           sce[[group]],
                                           statistics = c("sum", "prop.detected"),
                                           assay.type = "nRef",
                                           BPPARAM = BPPARAM)

    nalt <- scuttle::summarizeAssayByGroup(sce,
                                           sce[[group]],
                                           statistics = c("sum", "prop.detected"),
                                           assay.type = "nAlt",
                                           BPPARAM = BPPARAM)
   
    grp_pairs <- group_combinations(sce[[group]])
    
    stats <- BiocParallel::bplapply(seq_len(nrow(grp_pairs)),
                                    function(i) {
                                        gp <- grp_pairs[i, , drop = TRUE]
                                        ref <- assay(nref, "sum")[, c(gp$first, gp$second)]
                                        alt <- assay(nalt, "sum")[, c(gp$first, gp$second)]
                                        pvals <- calc_fisher_exact(ref, alt)
                                        ef <- alt / (ref + alt)
                                        d_editing_frequency <- ef[, 2] - ef[, 1]

                                        res <- data.frame(p.value = pvals,
                                                          dEF = d_editing_frequency,
                                                          row.names = rownames(ref))
                                        res
                                    }, BPPARAM = BPPARAM)
    res <- scran::combineMarkers(stats,
                                 grp_pairs,
                                 effect.field = "dEF",
                                 BPPARAM = BPPARAM,
                                 ...)

    depth_summary <- scran::summaryMarkerStats(sce,
                                               sce[[group]],
                                               assay.type = "depth",
                                               BPPARAM = BPPARAM)

    for(nm in names(res)) {
        de_stats <- res[[nm]]
        depths <- depth_summary[[nm]][rownames(de_stats), ]
        res[[nm]] <- cbind(depths, de_stats)
    }

    res
}

calc_fisher_exact <- function(ref, alt) {

    vals <- t(cbind(ref, alt))
    mode(vals) <- "integer"

    if(any(is.na(vals))) {
        cli::cli_abort("NA values not supported in fisher test")
    }

    stopifnot(nrow(vals) == 4)
    .Call(".fisher_exact", vals)
}

group_combinations <- function(x) {
    if(length(x) < 2) {
        cli::cli_abort("At least 2 groups must be present")
    }
    
    if(is.factor(x)) {
        grps <- levels(droplevels(x))
    } else {
        grps <- levels(as.factor(x))
    }
    
    grp_pairs <- expand.grid(grps, grps, stringsAsFactors = FALSE)
    colnames(grp_pairs) <- c("first", "second")
    grp_pairs <- grp_pairs[grp_pairs[, 1] != grp_pairs[, 2], ]
    grp_pairs <- grp_pairs[order(grp_pairs$first, grp_pairs$second), ]
    rownames(grp_pairs) <- NULL
    grp_pairs
}

