# Apply strand correction using gene annotations

Gene annotations are used to infer the likely strand of editing sites.
This function will operate on unstranded datasets which have been
processed using "unstranded" library type which reports variants with
respect to the + strand for all sites. The strand of the editing site
will be assigned the strand of overlapping features in the `genes_gr`
object. Sites with no-overlap, or overlapping features with conflicting
strands (+ and -) will be removed.

## Usage

``` r
correct_strand(rse, genes_gr)
```

## Arguments

- rse:

  RangedSummarizedExperiment object containing editing sites processed
  with "unstranded" setting

- genes_gr:

  GRanges object containing reference features to annotate the strand of
  the editing sites.

## Value

RangedSummarizedExperiment object containing pileup assays, with strand
corrected based on supplied genomic intervals.

## Examples

``` r
suppressPackageStartupMessages(library("GenomicRanges"))

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
fp <- FilterParam(library_type = "unstranded")
rse <- pileup_sites(bamfn, fafn, param = fp)

genes <- GRanges(c(
    "DHFR:200-400:+",
    "SPCS3:100-200:-",
    "SSR3:3-10:-",
    "SSR3:6-12:+"
))

correct_strand(rse, genes)
#> class: RangedSummarizedExperiment 
#> dim: 307 1 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(307): site_SSR3_3_1 site_SSR3_4_1 ... site_DHFR_399_1
#>   site_DHFR_400_1
#> rowData names(4): REF rpbz vdb sor
#> colnames(1): SRR5564269_Aligned.sortedByCoord.out.md.bam
#> colData names(1): sample
```
