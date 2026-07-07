# Filter out clustered sequence variants

Sequence variants of multiple allele types (e.g., `A -> G`, `A -> C`)
proximal to a putative editing site can be indicative of a region prone
to mis-alignment artifacts. Sites will be removed if variants of
multiple allele types are present within a given distance in genomic or
transcriptome coordinate space.

## Usage

``` r
filter_clustered_variants(
  rse,
  txdb,
  regions = c("transcript", "genome"),
  variant_dist = 100
)
```

## Arguments

- rse:

  [`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  containing editing sites

- txdb:

  [`GenomicFeatures::TxDb`](https://rdrr.io/pkg/GenomicFeatures/man/TxDb-class.html)

- regions:

  One of `transcript` or `genome`, specifying the coordinate system for
  calculating distances between variants.

- variant_dist:

  distance in nucleotides for determining clustered variants

## Value

[`SummarizedExperiment::SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with sites removed from object dependent on filtering applied.

## See also

Other se-filters:
[`filter_multiallelic()`](https://rnabioco.github.io/raer/reference/filter_multiallelic.md),
[`filter_splice_variants()`](https://rnabioco.github.io/raer/reference/filter_splice_variants.md)

## Examples

``` r
if(require("txdbmaker")){
  rse_adar_ifn <- mock_rse()
  rse <- rse_adar_ifn[seqnames(rse_adar_ifn) == "SPCS3"]

  # mock up a txdb with genes
  gr <- GRanges(c(
      "SPCS3:100-120:-",
      "SPCS3:325-350:-"
  ))
  gr$source <- "raer"
  gr$type <- "exon"
  gr$source <- NA
  gr$phase <- NA_integer_
  gr$gene_id <- c(1, 2)
  gr$transcript_id <- c("1.1", "2.1")
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)

  rse <- filter_multiallelic(rse)
  filter_clustered_variants(rse, txdb, variant_dist = 10)

}
#> Loading required package: txdbmaker
#> Loading required package: GenomicFeatures
#> Loading required package: AnnotationDbi
#> Warning: genome version information is not available for this TxDb object
#> ℹ `filter_multiallelic()`: removed 0 sites from 8 (8 remain)
#> ℹ `filter_clustered_variants()`:  removed 3 sites from 8 (5 remain)
#> class: RangedSummarizedExperiment 
#> dim: 5 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(5): site_SPCS3_99_1 site_SPCS3_227_1 site_SPCS3_241_1
#>   site_SPCS3_329_1 site_SPCS3_330_1
#> rowData names(5): REF rpbz vdb sor ALT
#> colnames(2): wt adar1_ko
#> colData names(1): sample
```
