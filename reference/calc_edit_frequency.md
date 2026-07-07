# Adds editing frequencies

Adds editing frequencies to an existing
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
object (created by
[`pileup_sites()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)).
The
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
with a new assay for editing frequencies for each site (`edit_freq`),
depth of coverage computed using the indicated edited nucleotides
(`depth`) and new `colData` columns with the number of edited sites
(`n_sites`) and the fraction of edits (`edit_idx`) is returned.

## Usage

``` r
calc_edit_frequency(
  rse,
  edit_from = "A",
  edit_to = "G",
  drop = FALSE,
  replace_na = TRUE,
  edit_frequency = 0,
  min_count = 1
)
```

## Arguments

- rse:

  A
  [RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object created by
  [`pileup_sites()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)

- edit_from:

  This should correspond to a nucleotide or assay (`A`, `C`, `G`, `T`,
  `Ref`, or `Alt`) you expect in the reference. Ex. for A to I editing
  events, this would be `A`.

- edit_to:

  This should correspond to a nucleotide or assay (`A`, `C`, `G`, `T`,
  `Ref`, or `Alt`) you expect in the editing site. Ex. for A to I
  editing events, this would be `G`.

- drop:

  If `TRUE`, the
  [RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  returned will only retain sites matching the specified `edit_from` and
  `edit_to` bases.

- replace_na:

  If `TRUE`, `NA` and `NaN` editing frequencies will be coerced to `0`.

- edit_frequency:

  The edit frequency cutoff used when calculating the number of sites.
  Set to `0` to require any non-zero editing frequency. The number of
  sites is stored as `n_sites` in the `colData`.

- min_count:

  The minimum number of reads required when enumerating number of
  editing sites detected.

## Value

[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
supplemented with `edit_freq` and `depth` assay.

## Examples

``` r
library(SummarizedExperiment)
rse_adar_ifn <- mock_rse()
rse <- calc_edit_frequency(rse_adar_ifn)
#> ℹ 6 sites had no coverage for calculating editing
assay(rse, "edit_freq")[1:5, ]
#>                        wt adar1_ko
#> site_SSR3_102_2 0.0000000        1
#> site_SSR3_125_2 1.0000000        0
#> site_SSR3_156_2 0.0000000        0
#> site_SSR3_176_2 0.6666667        0
#> site_SSR3_198_2 0.0400000        0
```
