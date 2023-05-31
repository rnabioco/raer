
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raer <a href="https://rnabioco.github.io/raer"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![R-CMD-check-bioc](https://github.com/rnabioco/raer/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/rnabioco/raer/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/raer/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/raer?branch=main)
<!-- badges: end -->

raer facilitates analysis of RNA adenosine editing in the
[Bioconductor](https://bioconductor.org/) ecosystem.

## Installation

You can install the development version of raer from
[GitHub](https://github.com/rnabioco/raer) with:

``` r
# install.packages("BiocManager")
BiocManager::install("rnabioco/raer")
```

## Quick start

raer provides methods to compute per site read count summaries from BAM
alignment files, either for known editing sites, or for all detected
sites.

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
#> rownames(1695): site_SSR3_1_2 site_SSR3_2_2 ... site_DHFR_517_2
#>   site_DHFR_518_2
#> rowData names(4): REF rpbz vdb sor
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
#>               ko wt
#> site_SSR3_1_2 13 12
#> site_SSR3_2_2 14 12
#> site_SSR3_3_2 14 12
#> site_SSR3_4_2 15 12
assays(rse)$nAlt[1:4, ]
#>               ko wt
#> site_SSR3_1_2  0  0
#> site_SSR3_2_2  0  0
#> site_SSR3_3_2  0  0
#> site_SSR3_4_2  0  0
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
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(2): ko wt
#> colData names(1): sample
```

`pileup_cells()` provides support for quantifying editing sites in
single cell libraries.

``` r
scbam_fn <- raer_example("5k_neuron_mouse_possort.bam")
outdir <- tempdir("sc_editing")

editing_sites <- GRanges(
    c(
        "2:579:-",
        "2:625:-",
        "2:589:-"
    ),
    REF = "A",
    ALT = "G"
)

cbs <- c(
    "CACCAAACAACAACAA-1",
    "TATTCCACACCCTCTA-1",
    "GACCTTCAGTTGTAAG-1"
)

sce <- pileup_cells(scbam_fn,
    sites = editing_sites,
    cell_barcodes = cbs,
    param = fp,
    output_directory = outdir
)
sce
#> class: SingleCellExperiment 
#> dim: 3 3 
#> metadata(0):
#> assays(2): nRef nAlt
#> rownames(3): site_2_579_2 site_2_625_2 site_2_589_2
#> rowData names(2): REF ALT
#> colnames(3): CACCAAACAACAACAA-1 TATTCCACACCCTCTA-1 GACCTTCAGTTGTAAG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

``` r
assays(sce)$nRef
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>              CACCAAACAACAACAA-1 TATTCCACACCCTCTA-1 GACCTTCAGTTGTAAG-1
#> site_2_579_2                  0                  0                  1
#> site_2_625_2                  0                  0                  0
#> site_2_589_2                  1                  1                  2
assays(sce)$nAlt
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>              CACCAAACAACAACAA-1 TATTCCACACCCTCTA-1 GACCTTCAGTTGTAAG-1
#> site_2_579_2                  1                  1                  1
#> site_2_625_2                  1                  1                  1
#> site_2_589_2                  0                  0                  0
```

## Related work

Core routines in `raer` are implemented using the `htslib` library and
methods from `samtools` and `bcftools`. `raer` builds off of approaches
from other rna editing detection tools:

- [REDItools](https://github.com/BioinfoUNIBA/REDItools) from [Picardi
  E, Pesole G](https://doi.org/10.1093/bioinformatics/btt287)  
- [JACUSA2](https://github.com/dieterich-lab/JACUSA2) from [Piechotta M
  et al](https://doi.org/10.1186/s12859-016-1432-8)  
- [deNovo-Detect](https://github.com/a2iEditing/deNovo-Detect) from
  [Gabey O et al](https://doi.org/10.1038/s41467-022-28841-4)  
- [RNAEditingIndexer](https://github.com/a2iEditing/RNAEditingIndexer)
  from [Roth SH et al](https://doi.org/10.1038/s41592-019-0610-9)  
- [SAILOR](https://github.com/YeoLab/sailor) from [Washburn MC et
  al](https://10.1016/j.celrep.2014.01.011)
