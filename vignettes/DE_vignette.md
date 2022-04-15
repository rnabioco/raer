# DE tutorial

Steps to perform differential editing using `Raer`

## Set up
First load in requried packages

```{r}
library(raer)
library(here)
library(tidyverse)
library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
```


Set directories and read in the data (eventually, we should make this downloadable)

```{r}
res_dir <- here("../working_dir/results/GSE99249/")

# Read in data
se <- readRDS(file.path(res_dir, "merged_se_with_meta.rds"))
```

## Preprocess
We start by using the `add_editing_frequencies` function to identify the percent of edits for each position and sample. Several pre-built options are available: A to I ("AI"), U to C ("UC"), C to U ("CU"), and G to A ("GA"). You can specify A to I editing by setting `type = AI` (this is also the default).

```{r}
se_filtered <- add_editing_frequencies(se, type = "AI",
                                       edit_frequency = 0.01,
                                       min_count = 10)
```

If you are interested in some other kind of editing, you can use a custom event. To do this, set `type = "none"` and `edit_from = nucleotide_1` and `edit_to = nucleotide_2`. Here, the `edit_from` should be the nucleotide in the reference (for A to I this would be `A`) and `edit_to` should be the nucleotide after editing (for A to I this would be `G`). Both `edit_from` and `edit_to` must be a valid nucleotide (A, T, C, G).

```{r}
se_filtered <- add_editing_frequencies(se, type = "none",
                                       edit_frequency = 0.01,
                                       min_count = 10,
                                       edit_from = "C", edit_to = "T")
```

You can also make plots that show the number of edited sites and the frequency of editing per sample. You can make these plots using `make_plots = TRUE`. If `make_plots = TRUE` and no save dir is specified, the plots and the filtered object will be returned as a list

```{r}
se_results <- add_editing_frequencies(se, type = "AI",
                                      edit_frequency = 0.01,
                                      min_count = 10,
                                      make_plots = TRUE)
                                      
se_filtered <- se_results[[1]]
site_plot <- se_results[[2]]
frequency_plot <- se_results[[3]]
```

If you would rather directly save these plots, you can provide a save directory. The plots will be saved as "number_of_sites.pdf" and "editing_index.pdf". If the plots are saved, only the filtered object will be returned.

```{r}
se_filtered <- add_editing_frequencies(se, type = "AI",
                                       edit_frequency = 0.01,
                                       min_count = 10,
                                       make_plots = TRUE,
                                       save_dir = "results")
```

Once your object has been filtered, you can prepare it for DE. This means making an assay of counts that contains both your alt and ref allele.

```{r}
deobj <- prep_for_de(se_filtered,
                     min_prop = 0.1,
                     max_prop = 0.9,
                     min_samples = 3)
```

## Run differential editing (DESeq2)

At this stage, you can use the object to perform DE yourself or you can continue with our pre built functions

For differential editing, we use the design `design <- ~0 + condition:sample + condition:count`.

For the samples, you can leave as is or combine so the same sample name shows up in both the treatment and control. These results are not identical but they are close. In my hands, the same genes come out, but the p values and log fold change values are slightly different.

It is probably best to update the levels of your object, but if you don't, this will still work.

To run using `DESEq2`, set `type = DESEq2` in the `de_results` function. This function requires you to specify what your control and treatment are from your condition column of your `deobj`.

```{r}
de_results <- perform_de(deobj, type = "DESeq2", sample_col = "sample",
                         condition_col = "genotype",
                         condition_control = "Wildtype",
                         condition_treatment = "ADAR1KO")
```

This returns a list of the dds object, the full results, the significant results, and the model matrix.


If you'd rather use the sample information as described in the `DESEq2` vignette, you can make your sample names go across your treatment and control

```{r}
sample_mapping = c("SRR5564260" = "s1",
                   "SRR5564270" = "s1",
                   "SRR5564261" = "s2",
                   "SRR5564271" = "s2",
                   "SRR5564268" = "s3",
                   "SRR5564274" = "s3",
                   "SRR5564269" = "s4",
                   "SRR5564275" = "s4",
                   "SRR5564272" = "s5",
                   "SRR5564276" = "s5",
                   "SRR5564273" = "s6",
                   "SRR5564277" = "s6")

deobj$sample_new <- sample_mapping[deobj$sample]

deobj$sample_new <- factor(deobj$sample_new)

de_results <- perform_de(deobj, type = "DESeq2", sample_col = "sample_new",
                         condition_col = "genotype",
                         condition_control = "Wildtype",
                         condition_treatment = "ADAR1KO")

```

## Run differential editing (edgeR)

At this stage, you can use the object to perform DE yourself or you can continue with our pre built functions

For differential editing, we use the design `design <- ~0 + condition:sample + condition:count`.

For the samples, you can leave as is or combine so the same sample name shows up in both the treatment and control. These results are not identical but they are close. In my hands, the same genes come out, but the p values and log fold change values are slightly different.

It is probably best to update the levels of your object, but if you don't, this will still work.

To run using `edgeR`, set `type = edgeR` in the `de_results` function. This function requires you to specify what your control and treatment are from your condition column of your `deobj`.

```{r}
de_edger_results <- perform_de(deobj, type = "edgeR",
                               sample_col = "sample",
                               condition_col = "genotype",
                               condition_control = "Wildtype",
                               condition_treatment = "ADAR1KO")
```

This returns a list of the dds object, the full results, the significant results, and the model matrix.


If you'd rather use the sample information as described in the `DESEq2` vignette, you can make your sample names go across your treatment and control

```{r}
sample_mapping = c("SRR5564260" = "s1",
                   "SRR5564270" = "s1",
                   "SRR5564261" = "s2",
                   "SRR5564271" = "s2",
                   "SRR5564268" = "s3",
                   "SRR5564274" = "s3",
                   "SRR5564269" = "s4",
                   "SRR5564275" = "s4",
                   "SRR5564272" = "s5",
                   "SRR5564276" = "s5",
                   "SRR5564273" = "s6",
                   "SRR5564277" = "s6")

deobj$sample_new <- sample_mapping[deobj$sample]

deobj$sample_new <- factor(deobj$sample_new)

de_edger_results <- perform_de(deobj, type = "edgeR",
                               sample_col = "sample_new",
                               condition_col = "genotype",
                               condition_control = "Wildtype",
                               condition_treatment = "ADAR1KO")

```


## Perform your own comparisons (DESeq2)

I have found it more straight forward to pull out comparisions of interest following [this tutorial](https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md). To ask for other comparisons with this dataset, you first need to identify the specific matrix design of the comparison.

The model matrix used for your experiment will be returned from the `perform_de` function. To get these results, first pull them out of the return list

```{r}
mod_mat <- de_results$model_matrix
```

We will also need to get the DESeq object from this list
```{r}
dds <- de_results$deseq_obj
```

We now can subset this matrix to the conditions we are interested in

```{r}
alt_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
                                  dds$count == "alt",])
ref_treatment <- colMeans(mod_mat[dds$condition == condition_treatment &
                                  dds$count == "ref",])

alt_control <- colMeans(mod_mat[dds$condition == condition_control &
                                dds$count == "alt",])
ref_control <- colMeans(mod_mat[dds$condition == condition_control &
                                dds$count == "ref",])
```

We can now use these to pull out specific comparisions using the `results` function.

We can look at the differences between treatment and control for just the alt allele:

```{r}
alt_treat_control <- results(dds, contrast = alt_treatment - alt_control)
```

Or just the ref allele:

```{r}
ref_treat_control <- results(dds, contrast = ref_treatment - ref_control)
```

Or just find the alt vs ref allele for treatment:
```{r}
alt_ref_treat <- results(dds, contrast = alt_treatment - ref_treatment)
```

To find the alt vs ref specific to the treatment (what we returned as the results from `perform_de`):

```{r}
  treatment_vs_control <- results(dds,
                                  contrast = (alt_treatment - ref_treatment) - 
                                    (alt_control - ref_control))
```

## Perform your own comparisons (edgeR)

I have found it more straight forward to pull out comparisions of interest following [this tutorial](https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md). To ask for other comparisons with this dataset, you first need to identify the specific matrix design of the comparison.

The model matrix used for your experiment will be returned from the `perform_de` function. To get these results, first pull them out of the return list

```{r}
mod_mat <- de_results$model_matrix
```

We will also need to get the DESeq object from this list
```{r}
deobj <- de_results$deseq_obj
```

We now can subset this matrix to the conditions we are interested in

```{r}
alt_treatment <- colMeans(mod_mat[deobj$condition == condition_treatment &
                                 deobj$count == "alt",])
ref_treatment <- colMeans(mod_mat[deobj$condition == condition_treatment &
                                 deobj$count == "ref",])

alt_control <- colMeans(mod_mat[deobj$condition == condition_control &
                               deobj$count == "alt",])
ref_control <- colMeans(mod_mat[deobj$condition == condition_control &
                               deobj$count == "ref",])
```

We can now use these to pull out specific comparisions using the `results` function.

We can look at the differences between treatment and control for just the alt allele:

```{r}
alt_treat_control <- glmLRT(deobj, contrast = alt_treatment - alt_control)

alt_treat_control <- topTags(alt_treat_control, n = nrow(alt_treat_control))
```

Or just the ref allele:

```{r}
ref_treat_control <- glmLRT(deobj, contrast = ref_treatment - ref_control)

ref_treat_control <- topTags(ref_treat_control, n = nrow(ref_treat_control))
```

Or just find the alt vs ref allele for treatment:
```{r}
alt_ref_treat <- glmLRT(deobj, contrast = alt_treatment - ref_treatment)

alt_ref_treat <- topTags(alt_ref_treat, n = nrow(alt_ref_treat))
```

To find the alt vs ref specific to the treatment (what we returned as the results from `perform_de`):

```{r}
treatment_vs_control <- glmLRT(deobj, contrast = (alt_treatment - ref_treatment) - 
                  (alt_control - ref_control))

treatment_vs_control <- topTags(treatment_vs_control,
                              n = nrow(treatment_vs_control))
```
