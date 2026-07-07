# Extract regions surrounding splice sites

Extract intervals at splice sites and their adjacent regions.

## Usage

``` r
get_splice_sites(txdb, slop = 4)
```

## Arguments

- txdb:

  [`GenomicFeatures::TxDb`](https://rdrr.io/pkg/GenomicFeatures/man/TxDb-class.html)

- slop:

  The number of bases upstream and downstream of splice site to extract

## Value

[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
containing positions of splice sites, with flanking bases.

## Examples

``` r
if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    res <- get_splice_sites(txdb)
    res[1:5]
}
#> Loading required package: TxDb.Hsapiens.UCSC.hg38.knownGene
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames      ranges strand
#>          <Rle>   <IRanges>  <Rle>
#>   [1]     chr1 11208-11215      +
#>   [2]     chr1 11208-11215      +
#>   [3]     chr1 11668-11675      +
#>   [4]     chr1 11668-11675      +
#>   [5]     chr1 11668-11675      +
#>   -------
#>   seqinfo: 711 sequences from an unspecified genome; no seqlengths
```
