
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

res <- get_pileup(bamfn, fafn, bedfile = bedfn)
res[1:5, ]
#> GRanges object with 5 ranges and 10 metadata columns:
#>       seqnames    ranges strand |         Ref         Var      nRef      nVar
#>          <Rle> <IRanges>  <Rle> | <character> <character> <integer> <integer>
#>   [1]     SSR3       201      - |           A           -        14         0
#>   [2]     SSR3       202      - |           T           -        14         0
#>   [3]     SSR3       203      - |           T           -        14         0
#>   [4]     SSR3       204      - |           T           -        15         0
#>   [5]     SSR3       205      - |           T           -        16         0
#>              nA        nT        nC        nG        nN        nX
#>       <integer> <integer> <integer> <integer> <integer> <integer>
#>   [1]        14         0         0         0         0         2
#>   [2]         0        14         0         0         0         2
#>   [3]         0        14         0         0         0         2
#>   [4]         0        15         0         0         0         1
#>   [5]         0        16         0         0         0         0
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```

Multiple BAM files can be processed, enabling rapid comparisons of
RNA-Seq vs. WGS or WXS data, or RNA-Seq vs RNA-seq (e.g., ADAR WT VS
ADAR KO).

The `FilterParam()` class holds multiple options for customizing the
output of `get_pileup()`.

To facilitate comparisons across groups, pileups can be combined in a
`RangedSummarizedExperiment`.

``` r
fp <- FilterParam(
  only_keep_variants = TRUE,
  library_type = "fr-first-strand",
  min_nucleotide_depth = 2
)

plps <- get_pileup(
  c(bam2fn, bamfn),
  fafn,
  filterParam = fp
)

create_se(plps)
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(7): Var nRef ... nC nG
#> rownames(74): SSR3_102_- SSR3_125_- ... DHFR_430_- DHFR_513_-
#> rowData names(1): Ref
#> colnames(2): sample_1 sample_2
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
