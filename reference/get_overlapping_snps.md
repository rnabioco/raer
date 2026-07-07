# Retrieve SNPs overlapping intervals

This function will find SNPs overlapping supplied intervals using a
SNPlocs package. The SNPs can be returned in memory (as GPos objects) or
written to disk as a bed-file (optionally compressed).

## Usage

``` r
get_overlapping_snps(gr, snpDb, output_file = NULL)
```

## Arguments

- gr:

  Intervals to query

- snpDb:

  A reference ot a SNPlocs database

- output_file:

  A path to an output file. If supplied the file can be optionally
  compressed by including a ".gz" suffix. If not supplied, SNPS will be
  returned as a
  [GenomicRanges::GPos](https://rdrr.io/pkg/GenomicRanges/man/GPos-class.html)
  object

## Value

GPos object containing SNPs overlapping supplied genomic intervals

## Examples

``` r
if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
    gr <- GRanges(rep("22", 10),
        IRanges(seq(10510077, 10610077, by = 1000)[1:10], width = 250),
        strand = "+"
    )
    get_overlapping_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
}
#> UnstitchedGPos object with 1 position and 0 metadata columns:
#>       seqnames       pos strand
#>          <Rle> <integer>  <Rle>
#>   [1]       22  10511116      *
#>   -------
#>   seqinfo: 25 sequences (1 circular) from GRCh38.p2 genome
```
