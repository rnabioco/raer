#' Adds editing frequencies
#'
#' Adds editing frequencies to an existing SummarizedExperimentobject (created by
#' `create_se`). The SummarizedExperiment with a new assay for editing frequences for each site
#'  (edit_freq) and new colData columns with the number of edited sites (n_sites) and the
#' fraction of edits (edit_idx) is returned.
#'
#' @param se_object A SummarizedExperiment object created by `create_se`
#' @param edit_from This should be a nucleotide (A, C, G, or T)
#' corresponding to the nucleotide you expect in the reference. Ex. for A to I
#' editing events, this would be "A". If NULL, then editing frequencies will be
#' calculated using the `nVar` and `nRef` values.
#' @param edit_to This should be a nucleotide (A, C, G, or T) and should
#' correspond to the nucleotide you expect after the editing event. Ex. for A to I
#' editing events, this would be "G". If NULL, then editing frequencies will be
#' calculated using the `nVar` and `nRef` values.
#' @param drop If TRUE, the summarizedExperiment returned will only retain sites
#' matching the specified `edit_from` and `edit_to` bases.
#' @param replace_na If TRUE, NA and NaN editing frequencies will be coerced to 0.
#' @param edit_frequency  The edit frequency cutoff used when calculating the
#' number of sites. Set to 0 to require any non-zero editing frequency. The number of sites is
#' stored as n_sites in the colData.
#' @param min_count The minimum number of reads required when enumerating number of
#' editing sites detected.
#'
#' @examples
#' example(create_se, echo = FALSE)
#' se <- calc_edit_frequency(se)
#' assay(se, "edit_freq")[1:5, ]
#'
#' @import SummarizedExperiment
#' @export
calc_edit_frequency <- function(se_object,
                                edit_from = NULL,
                                edit_to = NULL,
                                drop = FALSE,
                                replace_na = TRUE,
                                edit_frequency = 0,
                                min_count = 1){

  # Set edit to and from for pre defined types
  if (is.null(edit_from) | is.null(edit_to)){
    edit_from <- "Ref"
    edit_to <- "Var"
  } else if (!(edit_from %in% c("A", "C", "G", "T")) |
             !(edit_to %in% c("A", "C", "G", "T"))) {
    stop("`edit_to` and `edit_from` must be nucleotides!")
  }

  from_col <- paste0("n", edit_from)
  to_col <- paste0("n", edit_to)

  if(drop && from_col != "nRef"){
    se_object <- se_object[mcols(rowRanges(se_object))$Ref == edit_from, ]
  }

  assay(se_object, "edit_freq") <- assay(se_object, to_col) /
    (assay(se_object, from_col) + assay(se_object, to_col))

  if(replace_na){
    assay(se_object, "edit_freq")[is.na(assay(se_object, "edit_freq"))] <- 0
  }

  se_object <- count_edits(se_object, edit_frequency, min_count,
                           edit_from, edit_to)

  se_object
}

#' Counts edits
#'
#' Counts edits per sample and add new colData columns with the number of edited sites
#' (n_sites) and the  fraction of edits (edit_idx). This function should be called by
#' `calc_edit_frequency` and is not meant to be used directly.
#'
#' @param se_filtered A SummarizedExperiment object created by `create_se` and
#' processed by `calc_edit_frequency`
#' @param edit_from OPTIONAL if not using a pre-built type, you can specify
#' your own editing. This should be a nucleotide (A, C, G, or T) and should
#' correspond to the nucleotide you expect in the reference. Ex. for A to I
#' editing events, this would be "A". If type is not "AI", both edit from
#' and edit_to must be set.
#' @param edit_to OPTIONAL if not using a pre-built type, you can specify
#' your own editing. This should be a nucleotide (A, C, G, or T) and should
#' correspond to the nucleotide you expect after the editing event. Ex. for A to I
#' editing events, this would be "G". If type is not "AI", both edit from
#' and edit_to must be set.
#' @param edit_frequency OPTIONAL the edit frequency used to determine the number of sites.
#' Default is 0.01.
#' @param min_count OPTIONAL the number of reads used to determine the number of edited sites.
#' Default is 10.
#'
#' @import SummarizedExperiment
count_edits <- function(se_filtered, edit_frequency = 0.01, min_count = 10,
                         edit_from = NULL, edit_to = NULL){

  n_pass_filter <- colSums((assay(se_filtered, "edit_freq") > edit_frequency) &
                             ((assay(se_filtered, paste0("n", edit_from)) +
                                 assay(se_filtered, paste0("n", edit_to))) >=
                                min_count))

  colData(se_filtered)$n_sites <- n_pass_filter

  edit_idx <- colSums(assay(se_filtered, paste0("n", edit_to))) /
    (colSums(assay(se_filtered, paste0("n", edit_from))) +
       colSums(assay(se_filtered, paste0("n", edit_to))))

  colData(se_filtered)$edit_idx <- edit_idx

  return(se_filtered)
}


#' Makes summary plots of editing
#'
#' Generates plots of the number of sites edited per sample and the percent
#' of editing events per sample. This function is written to be called directly
#' by `calc_edit_frequency`
#'
#' @param se_object A SummarizedExperiment object created by `create_se`
#' @param colors OPTIONAL The colors of the replicates. If no colors are
#' provided, Set1 from `RColorBrewer will be used
#' @param meta_col The column in colData to be used to separate out samples
#' based on the condition. For example, "genotype", "treatment", or
#' "genotype_treatment." Default is "genotype_treatment."
#' @param replicate The column in colData contining information about the
#' replicates. Default is "rep".
#'
#' @import ggplot2
make_editing_plots <- function(se_object, colors = NULL,
                               meta_col = "genotype_treatment",
                               replicate = "rep"){
  if (is.null(colors)){
    if (!requireNamespace("RColorBrewer", quietly = TRUE)){
    stop(paste0("Package \"RColorBrewer\" needed for plotting if you don't provide colors.",
      " Please install it or provide colors with `colors = c()`."),
      call. = FALSE)
    }
    nColors <- length(unique(colData(se_object)[[replicate]]))
    if(nColors > 9){
      cols <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors)
    } else {
      cols <- RColorBrewer::brewer.pal(9, "Set1")
    }
  }
  plot_dat  <- as.data.frame(colData(se_object))
  names(plot_dat)[names(plot_dat) == meta_col] <- "sample_info"

  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes_string("sample_info", "n_sites")) +
    ggplot2::geom_col(ggplot2::aes(fill = rep),
                      position = ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL,
         y = "# of sites detected") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                              vjust = 0.5))

  p2 <- ggplot2::ggplot(plot_dat, ggplot2::aes_string("sample_info", "edit_idx")) +
    ggplot2::geom_col(ggplot2::aes(fill = rep),
                      position = ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL,
         y = "Editing Index") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  return(list(p1, p2))

}

#' Make summarized experiment object for DE

#' Generates a SummarizedExperiment object for use with edgeR or DESeq2
#' will generate a counts assay with a matrix formated with 2 columns per sample
#'
#' @param se A SummarizedExperiment object
#' @param type OPTIONAL the type of editing event to add. Currently, only
#' A to I is supported ("AI") which is the default, but your own custom can
#' be added by setting this to "none".
#' @param edit_from OPTIONAL if not using a pre-built type, you can specify
#' your own editing. This should be a nucleotide (A, C, G, or T) and should
#' correspond to the nucleotide you expect in the reference. Ex. for A to I
#' editing events, this would be "A". If type is not "AI", both edit from
#' and edit_to must be set.
#' @param edit_to OPTIONAL if not using a pre-built type, you can specify
#' your own editing. This should be a nucleotide (A, C, G, or T) and should
#' correspond to the nucleotide you expect after the editing event. Ex. for A to I
#' editing events, this would be "G". If type is not "AI", both edit from
#' and edit_to must be set.
#' @param min_prop OPTIONAL the min proporation of reads edited at a site.
#' At least min_samples need to pass this to keep the site. Default is 0.1.
#' @param max_prop OPTIONAL the max proporation of reads edited at a site.
#' At least min_samples need to pass this to keep the site. Default is 0.9.
#' @param min_samples OPTIONAL the minimum number of samples passing the cutoffs
#' to keep a site. Default is 3.
#'
#' @import SummarizedExperiment
#' @examples
#' example(create_se, echo = FALSE)
#' se <- calc_edit_frequency(se)
#' dse <- prep_for_de(se)
#' assay(dse, "counts")
#' dse
#' @export
prep_for_de <- function(se,
                        type = "AI",
                        edit_from = NULL, edit_to = NULL,
                        min_prop = 0.1,
                        max_prop = 0.9,
                        min_samples = 3){

  # Set edit to and from for pre defined types
  if(type == "AI"){
    edit_from <- "A"
    edit_to <- "G"
  } else if (is.null(edit_from) | is.null(edit_to)){
    stop("If not using a pre built type 'AI', `edit_from` and `edit_to` must be set")
  } else if (!(edit_from %in% c("A", "C", "G", "T")) | !(edit_to %in% c("A", "C", "G", "T"))) {
    stop("`edit_to` and `edit_from` must be nucleotides!")
  }

  # Only keep locations that pass cutoffs in a certain number of samples
  pass_cutoff <- (assay(se, "edit_freq") >= min_prop) &
    (assay(se, "edit_freq") <= max_prop)
  se <- se[rowSums(pass_cutoff) >= min_samples, ]

  # Set the ref and alternate allele and create a count table with both
  ref <- assay(se, paste0("n", edit_from))
  colnames(ref) <- paste0(colnames(ref), "_ref")
  alt <- assay(se, paste0("n", edit_to))
  colnames(alt) <- paste0(colnames(alt), "_alt")
  res <- cbind(ref, alt)
  mdata <- colData(se)

  # Join the meta data for all samples
  ref_mdata <- alt_mdata <- mdata
  rownames(ref_mdata) <- colnames(ref)
  ref_mdata$count <- "ref"
  rownames(alt_mdata) <- colnames(alt)
  alt_mdata$count <- "alt"
  mdata <- rbind(ref_mdata, alt_mdata)

  # create a new SummarizedExperiment
  res <- SummarizedExperiment(assays = list(counts = res),
                              colData = mdata)
  return(res)
}

#' Perform differential editing

#' Uses either edgeR or DESeq2 to perform differential editing analysis.
#' This will work for simple designs that have 1 treatment and 1 control.
#' For more complex designs, we suggest you perform your own.
#'
#' At the moment, this function will only find editing events specific
#' to the treatment, but it will be pretty straight forward to add other
#' possible return values.
#'
#' @param deobj A SummarizedExperiment object prepared for de by `prep_for_de`
#' @param type OPTIONAL if edgeR or DESeq should be run. Default is edgeR
#' @param sample_col OPTIONAL the name of the column from colData(deobj) that
#' contains your sample information. Default is sample. If you do not have a
#' column named "sample", you must provide the appropriate sample column
#' @param condition_col OPTIONAL the name of the column from colData(deobj) that
#' contains your treatment information. Default is condition, If you do not have a
#' column named "condition", you must provide the appropriate condition column
#' @param condition_control The name of the control condition. This must be a variable
#' in your condition_col of colData(deobj). No default provided.
#' @param condition_treatment The name of the treatment condition. This must be a variable
#' in your condition_col of colData(deobj).
#'
#' @examples
#' example(create_se, echo = FALSE)
#' se <- calc_edit_frequency(se)
#' dse <- prep_for_de(se)
#' res <- perform_de(dse, condition_control = "WT", condition_treatment = "KO")
#' res$sig_results[1:5, ]
#'
#' @returns A named list
#' - de_obj: The edgeR or deseq object used for differential editing analysis
#' - results_full: Unfiltered differenital editing results
#' - sig_results: Filtered differenial editing (FDR < 0.05)
#' - model_matrix: The model matrix used for generating DE results
#'
#' @import stringr
#' @export

perform_de <- function(deobj, type = "edgeR", sample_col = "sample",
                       condition_col = "condition",
                       condition_control = NULL,
                       condition_treatment = NULL){

  # Make sure all variables are present
  if (!sample_col %in% colnames(colData(deobj))){
    stop(paste0("somple_col must be a column in the colDat of your deobj. '",
                sample_col, "' not found in colnames(colData(deobj))!"))
  }
  if (!condition_col %in% colnames(colData(deobj))){
    stop(paste0("condition_col must be a column in the colData of your deobj. '",
                condition_col, "' not found in colnames(colData(deobj))!"))
  }

  # Rename columns based on the input
  if (sample_col != "sample"){
    colData(deobj)$sample <- NULL
  }
  if (condition_col != "condition"){
    colData(deobj)$condition <- NULL
  }

  new_columns  <- as.data.frame(colData(deobj))
  names(new_columns)[names(new_columns) == condition_col] <- "condition"
  names(new_columns)[names(new_columns) == sample_col] <- "sample"
  new_columns <- new_columns[c("sample", "condition")]

  colData(deobj) <- cbind(colData(deobj), new_columns)

  # Check that condition_control and condition_treatment are correct
  if (is.null(condition_control)){
    options <- unique(colData(deobj)$condition)
    stop(paste0("condition_control must be set. This should be the level of",
                " your meta data that corresponds to your control. Possible",
                " options from your experiment are: ",
                stringr::str_c(options, collapse = ", ")))
  }
  if (is.null(condition_treatment)){
    options <- unique(colData(deobj)$condition)
    stop(paste0("condition_treatment must be set. This should be the level of",
                " your meta data that corresponds to your control. Possible",
                " options from your experiment are: ",
                str_c(options, collapse = ", ")))
  }

  # Check that the treatment and control are in the object
  if (!condition_control %in% colData(deobj)$condition){
    options <- unique(colData(deobj)$condition)
    stop(paste0("condition_control must be a column in your deobj colData. '",
                condition_control, "' not found in the levels of the condition",
                " column of colData(deobj)! Possible",
                " options from your experiment are: ",
                str_c(options, collapse = ", ")))
  }
  if (!condition_treatment %in% colData(deobj)$condition){
    options <- unique(colData(deobj)$condition)
    stop(paste0("condition_treatment must be a column in your deobj colData. '",
                condition_treatment, "' not found in the levels of the condition",
                " column of colData(deobj)! Possible",
                " options from your experiment are: ",
                str_c(options, collapse = ", ")))
  }
  if (type == "edgeR"){
    results <- run_edger(deobj, condition_control, condition_treatment)
  } else if (type == "DESeq2"){
    results <- run_deseq2(deobj, condition_control, condition_treatment)
  } else {
    stop(paste0("Unrecognized type: '", type, "'. type must be either edgeR",
                " or DESeq2."))
  }
  return(results)
}

#' Perform differential editing with DESeq2

#' Uses DESeq2 to perform differential editing analysis.
#' This will work for simple designs that have 1 treatment and 1 control.
#' For more complex designs, we suggest you perform your own.
#' It will test if your sample column makes the model matrix not full
#' rank. If that happens, the model matrix will be modified to be full rank.
#' This is not intended to be called directly by the user, instead, this should
#' be called by `perform_de`
#'
#' At the moment, this function will only find editing events specific
#' to the treatment, but it will be pretty straight forward to add other
#' possible return values.
#'
#' @param deobj A SummarizedExperiment object prepared for de by `prep_for_de`
#' @param condition_control The name of the control condition. This must be a variable
#' in your condition_col of colData(deobj). No default provided.
#' @param condition_treatment The name of the treatment condition. This must be a variable
#' in your condition_col of colData(deobj).
#'
#' @importFrom stats model.matrix
run_deseq2 <- function(deobj, condition_control = NULL,
                       condition_treatment = NULL){
  if (!requireNamespace("DESeq2", quietly = TRUE)){
    stop(paste0("Package \"DESeq2\" needed to run differential analysis. Please install it."),
      call. = FALSE)
  }

  design <- ~0 + condition:sample + condition:count

  # See if the design is full rank, if not, remove sample info
  test_mat <- try(DESeq2::DESeqDataSetFromMatrix(countData = assay(deobj,
                                                                   "counts"),
                                         colData = colData(deobj),
                                         design = design), silent = TRUE)

  if (is(test_mat,"try-error")){
    sample <- deobj$sample
    condition <- deobj$condition
    count <- deobj$count

    design <- model.matrix(~0 + sample + condition:count)
    design <- design[,!grepl("countref", colnames(design))]
    mod_mat <- design
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = assay(deobj, "counts"),
                                colData = colData(deobj),
                                design = design)

  # We don't want size factors because we are looking at ratios within a sample
  DESeq2::sizeFactors(dds) <- rep(1, nrow(colData(deobj)))

  # TODO - figure out what modeling is best, also try local
  dds <- DESeq2::DESeq(dds, fitType = "parametric")

  if (!is(test_mat, "try-error")){
    mod_mat <- model.matrix(design(dds), colData(dds))
  }

  # Pull out the model matrix for all comparisons of interest
  alt_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
                                      dds$count == "alt",])
  ref_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
                                      dds$count == "ref",])

  alt_control <- colMeans(mod_mat[dds$condition == condition_control &
                                    dds$count == "alt",])
  ref_control <- colMeans(mod_mat[dds$condition == condition_control &
                                    dds$count == "ref",])

  # TODO add other possible return values and return as a list.
  # This finds editing specific to the condition
  treatment_vs_control <- DESeq2::results(dds,
                                          contrast = (alt_treatment -
                                                        ref_treatment) -
                                    (alt_control - ref_control))

  deseq_res <- as.data.frame(treatment_vs_control)
  deseq_res <- deseq_res[deseq_res$padj < 0.05, ]
  deseq_res <- deseq_res[order(deseq_res$log2FoldChange, decreasing = TRUE), ]

  return(list(de_obj = dds,
              results_full = treatment_vs_control,
              sig_results = deseq_res,
              model_matrix = mod_mat))
}

#' Perform differential editing with edgeR

#' Uses edgeR to perform differential editing analysis.
#' This will work for simple designs that have 1 treatment and 1 control.
#' For more complex designs, we suggest you perform your own.
#' It will test if your sample column makes the model matrix not full
#' rank. If that happens, the model matrix will be modified to be full rank.
#' This is not intended to be called directly by the user, instead, this should
#' be called by `perform_de`
#'
#' At the moment, this function will only find editing events specific
#' to the treatment, but it will be pretty straight forward to add other
#' possible return values.
#'
#' @param deobj A SummarizedExperiment object prepared for de by `prep_for_de`
#' @param condition_control The name of the control condition. This must be a variable
#' in your condition_col of colData(deobj). No default provided.
#' @param condition_treatment The name of the treatment condition. This must be a variable
#' in your condition_col of colData(deobj).
run_edger <- function(deobj, condition_control = NULL,
                       condition_treatment = NULL){

  if (!requireNamespace("edgeR", quietly = TRUE)){
    stop(paste0("Package \"edgeR\" needed to run differential analysis. Please install it."),
         call. = FALSE)
  }

  sample <- deobj$sample
  condition <- deobj$condition
  count <- deobj$count

  design <- model.matrix(~0 + condition:sample + condition:count)

  # Check if full rank, if not, fix
  if (!limma::is.fullrank(design)){
    design <- model.matrix(~0 + sample + condition:count)
    design <- design[,!grepl("countref", colnames(design))]
  }

  dge <- edgeR::DGEList(assay(deobj, "counts"), lib.size = rep(1, ncol(deobj)))
  dge <- edgeR::estimateDisp(dge)
  fit <- edgeR::glmFit(dge, design)

  # Pull out the model matrix for all comparisons of interest
  alt_treatment <- colMeans(design[deobj$condition == condition_treatment &
                                     deobj$count == "alt",])
  ref_treatment <- colMeans(design[deobj$condition == condition_treatment &
                                     deobj$count == "ref",])

  alt_control <- colMeans(design[deobj$condition == condition_control &
                                   deobj$count == "alt",])
  ref_control <- colMeans(design[deobj$condition == condition_control &
                                   deobj$count == "ref",])

  # TODO add other possible return values and return as a list.
  # This finds editing specific to the condition
  treatment_vs_control <- edgeR::glmLRT(fit, contrast = (alt_treatment - ref_treatment) -
                      (alt_control - ref_control))

  treatment_vs_control <- edgeR::topTags(treatment_vs_control,
                                  n = nrow(treatment_vs_control))

  edger_res <- as.data.frame(treatment_vs_control)
  edger_res <- edger_res[edger_res$FDR < 0.05, ]
  edger_res <- edger_res[order(edger_res$PValue), ]

  return(list(de_obj = fit,
              results_full = treatment_vs_control,
              sig_results = edger_res,
              model_matrix = design))

}
