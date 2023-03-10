
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
bam1fn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")

bams <- c("ko" = bam1fn, "wt" = bam2fn)

rse <- pileup_sites(bams, fafn)
```

To facilitate comparisons across groups, base count data and genomic
coordinates are stored in a `RangedSummarizedExperiment`.

``` r
suppressMessages(library(SummarizedExperiment))
rse
#> class: RangedSummarizedExperiment 
#> dim: 1695 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(1695): SSR3_1_- SSR3_2_- ... DHFR_517_- DHFR_518_-
#> rowData names(4): REF rbpz vpb sor
#> colnames(2): ko wt
#> colData names(1): sample
assays(rse)
#> List of length 7
#> names(7): ALT nRef nAlt nA nT nC nG
colData(rse)
#> DataFrame with 2 rows and 1 column
#>         sample
#>    <character>
#> ko          ko
#> wt          wt
```

``` r
assays(rse)$nRef[1:4, ]
#>          ko wt
#> SSR3_1_- 13 12
#> SSR3_2_- 14 12
#> SSR3_3_- 14 12
#> SSR3_4_- 15 12
assays(rse)$nAlt[1:4, ]
#>          ko wt
#> SSR3_1_-  0  0
#> SSR3_2_-  0  0
#> SSR3_3_-  0  0
#> SSR3_4_-  0  0
```

The `FilterParam()` class holds multiple options for customizing the
output of `pileup_sites()`.

``` r
fp <- FilterParam(
    only_keep_variants = TRUE,
    library_type = "fr-first-strand",
    min_depth = 2
)

rse <- pileup_sites(bams, fafn, param = fp)
rse
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(74): SSR3_102_- SSR3_125_- ... DHFR_430_- DHFR_513_-
#> rowData names(4): REF rbpz vpb sor
#> colnames(2): ko wt
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
