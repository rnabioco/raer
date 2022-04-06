#' Adds editing frequencies
#'
#' Adds editing frequencies to an existing SummarizedExperimentobject (created by
#' `create_se`). Currently A to I editing is supported as well as any custom editing
#' events. The SummarizedExperiment with a new assay for editing requences for each site
#'  (edit_freq) and new colData columns with the number of edited sites (n_sites) and the
#' fraction of edits (edit_idx) is returned. If make_plots is set to true and no save
#' directory is specified, a list including the SummarizedExperiment object and plots
#' of the number of sites editied and the fraction of editing events will be returned. 
#' 
#' @param se_object A SummarizedExperiment object created by `create_se`
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
#' @param make_plots OPTIONAL if plots showing the number and frequency of editing
#' events should be made. Default is False. If True, the colData columns containing
#' the treatment and replicate information should be provided. See `make_editing_plots`
#' @param save_dir OPTIONAL if the plots created should be saved, include a path to
#' the directory. The plots will be saved as "number_of_sites.pdf" and "editing_index.pdf".
#' If no directory is provided, the plots will be returned in a list with the
#' SummarizedExperiment object.
#' @param edit_frequency OPTIONAL the edit frequency used to determine the number of sites.
#' Default is 0.01.
#' @param min_count OPTIONAL the number of reads used to determine the number of edited sites.
#' Default is 10.
#' @param ... Options passed to `make_editing_plots`
#'
#' @import SummarizedExperiment
#' @export

# TODO look up other editing events
add_editing_frequencies <- function(se_object, type = "AI",
                                    edit_from = NULL, edit_to = NULL,
                                    make_plots = FALSE, save_dir = NULL,
                                    edit_frequency = 0.01, min_count = 10,
                                    ...){
  
  # Set edit to and from for pre defined types
  if(type == "AI"){
    edit_from <- "A"
    edit_to <- "G"
  } else if (is.null(edit_from) | is.null(edit_to)){
    stop("If not using a pre built type 'AI', `edit_from` and `edit_to` must be set")
  } else if (!(edit_from %in% c("A", "C", "G", "T")) | !(edit_to %in% c("A", "C", "G", "T"))) {
    stop("`edit_to` and `edit_from` must be nucleotides!")
  }
  
  # Only keep the sites with the edit of interest
  se_filtered <- se_object[mcols(rowRanges(se_object))$Ref == edit_from, ] # Use %in% if you want to try multiple
  assay(se_filtered, "edit_freq") <- assay(se_filtered, paste0("n", edit_to)) / 
    (assay(se_filtered, paste0("n", edit_from)) + 
       assay(se_filtered, paste0("n", edit_to)))
  assay(se_filtered, "edit_freq")[is.na(assay(se_filtered, "edit_freq"))] <- 0
  
  se_filtered <- count_edits(se_filtered, edit_frequency, min_count,
                              edit_from, edit_to)
  
  if(make_plots){
    plot_list <- make_editing_plots (se_filtered, ...)
    
    if(!is.null(save_dir)){
      cowplot::save_plot(file.path(save_dir, "number_of_sites.pdf"),
                         plot_list[[1]], base_asp = 1)
      
      cowplot::save_plot(file.path(save_dir, "editing_index.pdf"),
                         plot_list[[2]], base_asp = 1)
      
      return(se_filtered)
    } else {
      return(c(se_filtered, plot_list))
    }
  } else {
    return(se_filtered)
  }
}

#' Counts edits
#'
#' Counts edits per sample and add new colData columns with the number of edited sites
#' (n_sites) and the  fraction of edits (edit_idx). This function should be called by
#' `add_editing_frequencies` and is not meant to be used directly.
#' 
#' @param se_filtered A SummarizedExperiment object created by `create_se` and
#' processed by `add_editing_frequencies`
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
#' @export

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
#' by `add_editing_frequencies`
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
#' @import SummarizedExperiment
#' @export

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
  p1 <- colData(se_object) %>%
    as.data.frame() %>%
    dplyr::rename(sample_info = dplyr::all_of(meta_col)) %>%
    ggplot2::ggplot(ggplot2::aes(sample_info, n_sites)) +
    ggplot2::geom_col(ggplot2::aes(fill = rep),
                      position = ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL,
         y = "# of sites detected") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                              vjust = 0.5))
  
  p2 <- colData(se_object) %>%
    as.data.frame() %>%
    dplyr::rename(sample_info = dplyr::all_of(meta_col)) %>%
    ggplot2::ggplot(ggplot2::aes(sample_info, edit_idx)) +
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
#' @param sample_column OPTIONAL the name of the column from colData(deobj) that
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
#' @import SummarizedExperiment
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
  new_columns <- colData(deobj) %>%
    data.frame %>% 
    dplyr::rename(condition = condition_col) %>%
    dplyr::rename(sample = sample_col) %>%
    dplyr::select(c(sample, condition))
  
  colData(deobj) <- cbind(colData(deobj), new_columns)
  
  # Check that condition_control and condition_treatment are correct
  if (is.null(condition_control)){
    options <- unique(colData(deobj)$condition)
    stop(paste0("condition_control must be set. This should be the level of",
                " your meta data that corresponds to your control. Possible",
                " options from your experiment are: ",
                str_c(options, collapse = ", ")))
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
#' @param type OPTIONAL if EdgeR or DESeq should be run. Default is EdgeR
#' @param condition_control The name of the control condition. This must be a variable
#' in your condition_col of colData(deobj). No default provided.
#' @param condition_treatment The name of the treatment condition. This must be a variable
#' in your condition_col of colData(deobj).
#'
#' @import SummarizedExperiment
#' @export

run_deseq2 <- function(deobj, condition_control = NULL,
                       condition_treatment = NULL){
  if (!requireNamespace("DESeq2", quietly = TRUE)){
    stop(paste0("Package \"DESeq2\" needed to run differential analysis. Please install it."),
      call. = FALSE)
  }
  
  design <- ~0 + condition:sample + condition:count
  
  # See if the design is full rank, if not, remove sample info
  test_mat <- try(DESeqDataSetFromMatrix(countData = assay(deobj, "counts"),
                                         colData = colData(deobj),
                                         design = design), silent = TRUE)
  
  if (class(test_mat) == "try-error"){
    sample <- deobj$sample
    condition <- deobj$condition
    count <- deobj$count
    
    design <- model.matrix(~0 + sample + condition:count)
    design <- design[,!grepl("countref", colnames(design))]
    mod_mat <- design
  }
  
  dds <- DESeqDataSetFromMatrix(countData = assay(deobj, "counts"),
                                colData = colData(deobj),
                                design = design)
  
  # We don't want size factors because we are looking at ratios within a sample
  sizeFactors(dds) <- rep(1, nrow(colData(deobj)))
  
  # TODO - figure out what modeling is best, also try local
  dds <- DESeq(dds, fitType = "parametric")
  
  if (class(test_mat) != "try-error"){
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
  treatment_vs_control <- results(dds,
                                  contrast = (alt_treatment - ref_treatment) - 
                                    (alt_control - ref_control))
  
  deseq_res <- treatment_vs_control %>%
    data.frame %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
  
  return(list(deseq_obj = dds,
              results_full = treatment_vs_control,
              sig_results = deseq_res))
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
#' @param type OPTIONAL if EdgeR or DESeq should be run. Default is EdgeR
#' @param condition_control The name of the control condition. This must be a variable
#' in your condition_col of colData(deobj). No default provided.
#' @param condition_treatment The name of the treatment condition. This must be a variable
#' in your condition_col of colData(deobj).
#'
#' @import SummarizedExperiment
#' @export

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
  
  dge <- DGEList(assay(deobj, "counts"), lib.size = rep(1, ncol(deobj)))
  dge <- estimateDisp(dge)
  fit <- glmFit(dge, design)
  
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
  treatment_vs_control <- glmLRT(fit, contrast = (alt_treatment - ref_treatment) - 
                      (alt_control - ref_control))
  
  treatment_vs_control <- topTags(treatment_vs_control,
                                  n = nrow(treatment_vs_control))
  
  edger_res <- treatment_vs_control %>%
    data.frame %>%
    dplyr::filter(FDR < 0.05) %>%
    dplyr::arrange(PValue)
  
  return(list(deseq_obj = fit,
              results_full = treatment_vs_control,
              sig_results = edger_res))
  
}
