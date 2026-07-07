# Calculate confidence score for observing editing

Calculate a confidence score based on a Bayesian inverse probability
model as described by Washburn et al. Cell Reports. 2015, and
implemented in the SAILOR pipeline.

## Usage

``` r
calc_confidence(
  se,
  edit_to = "G",
  edit_from = "A",
  per_sample = FALSE,
  exp_fraction = 0.01,
  alpha = 0L,
  beta = 0L
)
```

## Arguments

- se:

  [`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  containing editing sites

- edit_to:

  edited base

- edit_from:

  non-edited base

- per_sample:

  if TRUE, calculate confidence per sample, otherwise edited and
  non-edited counts will be summed across all samples.

- exp_fraction:

  Numeric value between 0 and 1, specifying the expected error rate

- alpha:

  Pseudo-count to add to non-edited base counts

- beta:

  Pseudo-count to add to edited base counts

## Value

[`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with either a new assay or rowData column named "confidence" depending
on whether confidence is calculated `per_sample`.

## References

Washburn MC, Kakaradov B, Sundararaman B, Wheeler E, Hoon S, Yeo GW,
Hundley HA. The dsRBP and inactive editor ADR-1 utilizes dsRNA binding
to regulate A-to-I RNA editing across the C. elegans transcriptome. Cell
Rep. 2014 Feb 27;6(4):599-607. doi: 10.1016/j.celrep.2014.01.011. Epub
2014 Feb 6. PMID: 24508457; PMCID: PMC3959997.

SAILOR pipeline: https://github.com/YeoLab/sailor

## Examples

``` r
rse_adar_ifn <- mock_rse()
calc_confidence(rse_adar_ifn)
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(5): REF rpbz vdb sor confidence
#> colnames(2): wt adar1_ko
#> colData names(1): sample
calc_confidence(rse_adar_ifn, per_sample = TRUE)
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(8): ALT nRef ... nG confidence
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(2): wt adar1_ko
#> colData names(1): sample
```
