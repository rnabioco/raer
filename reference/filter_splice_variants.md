# Filter out sites near splice sites

Remove editing sites found in regions proximal to annotated splice
junctions.

## Usage

``` r
filter_splice_variants(rse, txdb, splice_site_dist = 4, ignore.strand = FALSE)
```

## Arguments

- rse:

  [`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with editing sites

- txdb:

  [`GenomicFeatures::TxDb`](https://rdrr.io/pkg/GenomicFeatures/man/TxDb-class.html)

- splice_site_dist:

  distance to splice site

- ignore.strand:

  if `TRUE`, ignore strand when comparing editing sites to splice sites

## Value

[`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with sites adjacent to splice sites removed.

## See also

Other se-filters:
[`filter_clustered_variants()`](https://rnabioco.github.io/raer/reference/filter_clustered_variants.md),
[`filter_multiallelic()`](https://rnabioco.github.io/raer/reference/filter_multiallelic.md)

## Examples

``` r
if(require("txdbmaker")) {
  rse_adar_ifn <- mock_rse()

  # mock up a txdb with genes
  gr <- GRanges(c(
      "DHFR:310-330:-",
      "DHFR:410-415:-",
      "SSR3:100-155:-",
      "SSR3:180-190:-"
  ))
  gr$source <- "raer"
  gr$type <- "exon"
  gr$source <- NA
  gr$phase <- NA_integer_
  gr$gene_id <- c(1, 1, 2, 2)
  gr$transcript_id <- rep(c("1.1", "2.1"), each = 2)
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)

  filter_splice_variants(rse_adar_ifn, txdb)
}
#> Warning: genome version information is not available for this TxDb object
#> ℹ `filter_splice_variants()`:  removed 5 sites from 74 (69 remain)
#> class: RangedSummarizedExperiment 
#> dim: 69 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(69): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(2): wt adar1_ko
#> colData names(1): sample

```
