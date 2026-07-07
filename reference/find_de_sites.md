# Perform differential editing

Use `edgeR` or `DESeq2` to perform differential editing analysis. This
will work for designs that have 1 treatment and 1 control group. For
more complex designs, we suggest you perform your own modeling.

## Usage

``` r
find_de_sites(
  deobj,
  test = c("edgeR", "DESeq2"),
  sample_col = "sample",
  condition_col = "condition",
  condition_control = NULL,
  condition_treatment = NULL
)
```

## Arguments

- deobj:

  A
  [RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object prepared for differential editing analysis by
  [`make_de_object()`](https://rnabioco.github.io/raer/reference/make_de_object.md)

- test:

  Indicate if `edgeR` or `DESeq2` should be run.

- sample_col:

  The name of the column from `colData(deobj)` that contains your sample
  information. Default is sample. If you do not have a column named
  "sample", you must provide the appropriate sample column

- condition_col:

  The name of the column from `colData(deobj)` that contains your
  treatment information. Default is condition, If you do not have a
  column named "condition", you must provide the appropriate condition
  column

- condition_control:

  The name of the control condition. This must be a variable in your
  `condition_col` of `colData(deobj)`. No default provided.

- condition_treatment:

  The name of the treatment condition. This must be a variable in your
  `condition_col` of `colData(deobj)`.

## Value

A named list:

- `de_obj`: The `edgeR` or `deseq` object used for differential editing
  analysis

- `results_full`: Unfiltered differential editing results

- `sig_results`: Filtered differential editing (FDR \< 0.05)

- `model_matrix`: The model matrix used for generating DE results

## Examples

``` r
library(SummarizedExperiment)
bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")

bams <- rep(c(bamfn, bam2fn), each = 3)
sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
names(bams) <- sample_ids

fp <- FilterParam(only_keep_variants = TRUE)
rse <- pileup_sites(bams, fafn, param = fp)
rse$condition <- substr(rse$sample, 1, 2)

rse <- calc_edit_frequency(rse)
#> ℹ 6 sites had no coverage for calculating editing
dse <- make_de_object(rse)
res <- find_de_sites(dse,
    condition_control = "WT",
    condition_treatment = "KO"
)
#> Using classic mode.
res$sig_results[1:3, ]
#>                      logFC   logCPM        LR       PValue        FDR
#> site_SSR3_244_2 -11.108016 21.06907 12.789217 0.0003486230 0.02314789
#> site_DHFR_361_2 -11.147610 21.96719 11.540891 0.0006808202 0.02314789
#> site_SSR3_254_2  -7.938091 21.10149  9.672241 0.0018707300 0.04240321
```
