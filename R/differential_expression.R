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

