# Generate a small RangedSummarizedExperiment object for tests and examples

A RangedSummarizedExperiment containing a subset of data from an RNA-seq
experiment to measure the effects of IFN treatment of cell lines with
wild-type or ADAR1-KO.

## Usage

``` r
mock_rse()
```

## Source

<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386593>

## Value

RangedSummarizedExperiment populated with pileup data

## References

<https://pubmed.ncbi.nlm.nih.gov/29395325/>

## Examples

``` r
mock_rse()
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(2): wt adar1_ko
#> colData names(1): sample
```
