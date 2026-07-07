# Filter out multi-allelic sites

Remove sites with multiple variant bases from a `SummarizedExperiment`.
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
gains a new column, `ALT`, that contains the variant allele detected at
each site.

## Usage

``` r
filter_multiallelic(se)
```

## Arguments

- se:

  [`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

## Value

[`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with multiallelic sites removed. A new column,`ALT` will be added to
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
indicating the single allele present at the site.

## See also

Other se-filters:
[`filter_clustered_variants()`](https://rnabioco.github.io/raer/reference/filter_clustered_variants.md),
[`filter_splice_variants()`](https://rnabioco.github.io/raer/reference/filter_splice_variants.md)

## Examples

``` r
rse_adar_ifn <- mock_rse()
filter_multiallelic(rse_adar_ifn)
#> ℹ `filter_multiallelic()`: removed 2 sites from 74 (72 remain)
#> class: RangedSummarizedExperiment 
#> dim: 72 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(72): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(5): REF rpbz vdb sor ALT
#> colnames(2): wt adar1_ko
#> colData names(1): sample
```
