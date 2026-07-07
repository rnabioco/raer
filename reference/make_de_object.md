# Make summarized experiment object for differential editing analysis

Generates a
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
object for use with `edgeR` or `DESeq2` . Will generate a `counts` assay
with a matrix formatted with 2 columns per sample, representing the
reference and editing allele counts.

## Usage

``` r
make_de_object(
  rse,
  edit_from = "A",
  edit_to = "G",
  min_prop = 0,
  max_prop = 1,
  min_samples = 1
)
```

## Arguments

- rse:

  A
  [RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object

- edit_from:

  This should correspond to a nucleotide or assay (`A`, `C`, `G`, `T`,
  `Ref`, or `Alt`) you expect in the reference. Ex. for A to I editing
  events, this would be `A`.

- edit_to:

  This should correspond to a nucleotide or assay (`A`, `C`, `G`, `T`,
  `Ref`, or `Alt`) you expect in the editing site. Ex. for A to I
  editing events, this would be `G`.

- min_prop:

  The minimum required proportion of reads edited at a site. At least
  `min_samples` need to pass this to keep the site.

- max_prop:

  The maximum allowable proportion of reads edited at a site. At least
  `min_samples` need to pass this to keep the site.

- min_samples:

  The minimum number of samples passing the `min_prop` and `max_prop`
  cutoffs to keep a site.

## Value

[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
for use with `edgeR` or `DESeq2`. Contains a `counts` assay with a
matrix formatted with 2 columns per sample (ref and alt counts).

## Examples

``` r
library(SummarizedExperiment)
rse_adar_ifn <- mock_rse()
rse <- calc_edit_frequency(rse_adar_ifn)
#> ℹ 6 sites had no coverage for calculating editing
dse <- make_de_object(rse, min_samples = 1)
assay(dse, "counts")[1:5, ]
#>                 wt_ref adar1_ko_ref wt_alt adar1_ko_alt
#> site_SSR3_102_2      0            0      0            1
#> site_SSR3_125_2      0            0      1            0
#> site_SSR3_156_2      1            0      0            0
#> site_SSR3_176_2      8           16     16            0
#> site_SSR3_198_2     24           15      1            0
dse
#> class: RangedSummarizedExperiment 
#> dim: 68 4 
#> metadata(0):
#> assays(1): counts
#> rownames(68): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_399_2
#>   site_DHFR_430_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(4): wt_ref adar1_ko_ref wt_alt adar1_ko_alt
#> colData names(4): sample n_sites edit_idx count
```
