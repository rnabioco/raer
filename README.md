
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raer <a href="https://rnabioco.github.io/raer"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![R-CMD-check-bioc](https://github.com/rnabioco/raer/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/rnabioco/raer/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/raer/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/raer?branch=main)
<!-- badges: end -->

raer facilitates analysis of RNA adenosine editing in the
[Bioconductor](https://bioconductor.org/) ecosystem.

🚧 **raer is under active development and functionality may change** 🚧

## Installation

You can install the development version of raer from
[GitHub](https://github.com/rnabioco/raer) with:

``` r
# install.packages("BiocManager")
BiocManager::install("rnabioco/raer")
```

## Quick start

raer provides methods to compute per site read count summaries from BAM
alignment files, either for known sites, or for all detected sites.

``` r
library(raer)
bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
bedfn <- raer_example("regions.bed")

rse <- get_pileup(c(bam2fn, bamfn), fafn, bedfile = bedfn)
```

To facilitate comparisons across groups, base count data and genomic
coordinates are stored in a `RangedSummarizedExperiment`.

``` r
suppressMessages(library(SummarizedExperiment))
rse
#> class: RangedSummarizedExperiment 
#> dim: 182 2 
#> metadata(0):
#> assays(7): Var nRef ... nC nG
#> rownames(182): SSR3_201_- SSR3_202_- ... DHFR_17_- DHFR_18_-
#> rowData names(3): Ref rbpz vpb
#> colnames(2): SRR5564277_Aligned.sortedByCoord.out.md.bam
#>   SRR5564269_Aligned.sortedByCoord.out.md.bam
#> colData names(1): sample
assays(rse)
#> List of length 7
#> names(7): Var nRef nVar nA nT nC nG
colData(rse)
#> DataFrame with 2 rows and 1 column
#>                                                             sample
#>                                                        <character>
#> SRR5564277_Aligned.sortedByCoord.out.md.bam SRR5564277_Aligned.s..
#> SRR5564269_Aligned.sortedByCoord.out.md.bam SRR5564269_Aligned.s..
rowRanges(rse)
#> GRanges object with 182 ranges and 3 metadata columns:
#>               seqnames    ranges strand |         Ref      rbpz       vpb
#>                  <Rle> <IRanges>  <Rle> | <character> <numeric> <numeric>
#>    SSR3_201_-     SSR3       201      - |           A -0.515758  0.036923
#>    SSR3_202_-     SSR3       202      - |           T  0.000000  0.000000
#>    SSR3_203_-     SSR3       203      - |           T  0.000000  0.000000
#>    SSR3_204_-     SSR3       204      - |           T  0.000000  0.000000
#>    SSR3_205_-     SSR3       205      - |           T  0.000000  0.000000
#>           ...      ...       ...    ... .         ...       ...       ...
#>   SPCS3_299_+    SPCS3       299      + |           C         0         0
#>   SPCS3_300_+    SPCS3       300      + |           T         0         0
#>     DHFR_16_-     DHFR        16      - |           G         0         0
#>     DHFR_17_-     DHFR        17      - |           G         0         0
#>     DHFR_18_-     DHFR        18      - |           A         0         0
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```

The `FilterParam()` class holds multiple options for customizing the
output of `get_pileup()`.

``` r
fp <- FilterParam(
  only_keep_variants = TRUE,
  library_type = "fr-first-strand",
  min_nucleotide_depth = 2
)

res <- get_pileup(bamfn, fafn, filterParam = fp)
res
#> class: RangedSummarizedExperiment 
#> dim: 31 1 
#> metadata(0):
#> assays(7): Var nRef ... nC nG
#> rownames(31): SSR3_102_- SSR3_228_- ... DHFR_430_- DHFR_513_-
#> rowData names(3): Ref rbpz vpb
#> colnames(1): SRR5564269_Aligned.sortedByCoord.out.md.bam
#> colData names(1): sample
```

## Related work

raer functionality builds off of previous work:

- Python package: [REDItools](https://github.com/BioinfoUNIBA/REDItools)
  from [Picardi E, Pesole
  G](https://doi.org/10.1093/bioinformatics/btt287)  
- Java tool: [JACUSA2](https://github.com/dieterich-lab/JACUSA2) from
  [Piechotta M et al](https://doi.org/10.1186/s12859-016-1432-8)  
- Python-based pipeline:
  [deNovo-Detect](https://github.com/a2iEditing/deNovo-Detect) from
  [Gabey O et al](https://doi.org/10.1038/s41467-022-28841-4)  
- Java-based tool:
  [RNAEditingIndexer](https://github.com/a2iEditing/RNAEditingIndexer)
  from [Roth SH et al](https://doi.org/10.1038/s41592-019-0610-9)
