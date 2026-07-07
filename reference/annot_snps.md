# Annotate known SNP positions

This function will annotate a
[GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) or
the rowRanges of a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with SNPs from a SNP package.

## Usage

``` r
annot_snps(obj, ...)

# S3 method for class 'GRanges'
annot_snps(
  obj,
  dbsnp,
  chrom = NULL,
  col_to_aggr = "RefSNP_id",
  drop = FALSE,
  genome = NULL,
  RLE = TRUE,
  ...
)

# S3 method for class 'SummarizedExperiment'
annot_snps(obj, ...)
```

## Arguments

- obj:

  GRanges or SummarizedExperiment object

- ...:

  For the generic, further arguments to pass to specific methods. Unused
  for now.

- dbsnp:

  SNPlocs package, see available packages from
  [`BSgenome::available.SNPs()`](https://rdrr.io/pkg/BSgenome/man/injectSNPs.html)

- chrom:

  only operate on a specified chromosome

- col_to_aggr:

  column from SNPlocs package to add to input. If multiple SNPs overlap
  these values will be concatenated as comma separated values.

- drop:

  If TRUE, remove sites overlapping SNPs

- genome:

  A BSgenome object, which if supplied, will be used to provide
  additional `snp_ref_allele` and `snp_alt_alleles` columns containing
  the reference and alt allele sequences, with respect to the positive
  strand. Additionally the snp sequences will be checked against the
  allele at the site if a column named `ALT` is present in object. The
  strand of the site will be used to determine if the `ALT` allele needs
  to be complemented prior to comparing against the SNP db (which always
  returns sequences w.r.t the plus strand).

- RLE:

  If TRUE, columns added will returned as
  [`S4Vectors::Rle()`](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html)
  vectors to reduce memory usage.

## Value

Either a GRanges or SummarizedExperiment object with a new column added
with information from `col_to_aggr` and optionally `snp_ref_allele`,
`snp_alt_alleles`, and `snp_matches_site` annotations.

## See also

[SNPlocs.Hsapiens.dbSNP144.GRCh38](https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh38.html)

## Examples

``` r
if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
    gr <- GRanges(rep("22", 10),
        IRanges(
            seq(10510077,
                10610077,
                by = 1000
            )[1:10],
            width = 250
        ),
        strand = "+"
    )
    genome(gr) <- "GRCh38.p2"
    annot_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
}
#> Loading required package: SNPlocs.Hsapiens.dbSNP144.GRCh38
#> Loading required package: BSgenome
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
#> Loading required package: BiocIO
#> Loading required package: rtracklayer
#> Warning: replacing previous import ‘utils::data’ by ‘BiocGenerics::data’ when loading ‘SNPlocs.Hsapiens.dbSNP144.GRCh38’
#> Warning: replacing previous import ‘utils::findMatches’ by ‘S4Vectors::findMatches’ when loading ‘SNPlocs.Hsapiens.dbSNP144.GRCh38’
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames            ranges strand | RefSNP_id
#>           <Rle>         <IRanges>  <Rle> |     <Rle>
#>    [1]       22 10510077-10510326      + |          
#>    [2]       22 10511077-10511326      + | rs4022986
#>    [3]       22 10512077-10512326      + |          
#>    [4]       22 10513077-10513326      + |          
#>    [5]       22 10514077-10514326      + |          
#>    [6]       22 10515077-10515326      + |          
#>    [7]       22 10516077-10516326      + |          
#>    [8]       22 10517077-10517326      + |          
#>    [9]       22 10518077-10518326      + |          
#>   [10]       22 10519077-10519326      + |          
#>   -------
#>   seqinfo: 1 sequence from GRCh38.p2 genome; no seqlengths
```
