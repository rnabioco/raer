
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raer

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/raer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rnabioco/raer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

raer is an R package that facilitates rapid interactive analysis of RNA
editing in R in the bioconductor ecosystem.

**raer is under active development and functionality may change**

## Installation

You can install the development version of raer from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rnabioco/raer")
```

## Quick start

The raer package provides methods to compute per site read count
summaries from bam files, either for known sites, or for all detected
sites.

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

The `FilterParam()` class holds multiple options for customizing the
output of `get_pileup()`.

``` r
fp <- FilterParam(only_keep_variants = TRUE)
res <- get_pileup(bamfn, fafn, filterParam = fp)
res
#> GRanges object with 31 ranges and 10 metadata columns:
#>        seqnames    ranges strand |         Ref         Var      nRef      nVar
#>           <Rle> <IRanges>  <Rle> | <character> <character> <integer> <integer>
#>    [1]     SSR3       102      - |           T          TG        12         1
#>    [2]     SSR3       228      - |           A          AG        13         1
#>    [3]     SSR3       244      - |           A          AG        18         1
#>    [4]     SSR3       254      - |           A          AG        18         1
#>    [5]     SSR3       258      - |           G          GA         8        10
#>    ...      ...       ...    ... .         ...         ...       ...       ...
#>   [27]     DHFR       300      - |           A          AG        56         1
#>   [28]     DHFR       332      - |           A          AG        55         1
#>   [29]     DHFR       336      - |           G          GT        51         1
#>   [30]     DHFR       430      - |           A          AT        37         1
#>   [31]     DHFR       513      - |           T          TC        35         1
#>               nA        nT        nC        nG        nN        nX
#>        <integer> <integer> <integer> <integer> <integer> <integer>
#>    [1]         0        12         0         1         0         6
#>    [2]        13         0         0         1         0         4
#>    [3]        18         0         0         1         0         1
#>    [4]        18         0         0         1         0         1
#>    [5]        10         0         0         8         0         2
#>    ...       ...       ...       ...       ...       ...       ...
#>   [27]        56         0         0         1         0         4
#>   [28]        55         0         0         1         0         2
#>   [29]         0         1         0        51         0         5
#>   [30]        37         1         0         0         0         5
#>   [31]         0        35         1         0         0         0
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```

Multiple bam files can be processed, which enables rapid comparisons of
RNA-Seq vs. WGS or WXS data, or RNA-Seq vs RNA-seq (ADAR WT VS ADAR KO).

``` r
fp <- FilterParam(
  only_keep_variants = TRUE,
  library_type = "fr-first-strand",
  min_nucleotide_depth = 2
)

plps <- get_pileup(c(bam2fn, bamfn),
  fafn,
  filterParam = fp
)
plps
#> [[1]]
#> GRanges object with 74 ranges and 10 metadata columns:
#>        seqnames    ranges strand |         Ref         Var      nRef      nVar
#>           <Rle> <IRanges>  <Rle> | <character> <character> <integer> <integer>
#>    [1]     SSR3       102      - |           T           -        15         0
#>    [2]     SSR3       125      - |           C          CG        21         1
#>    [3]     SSR3       156      - |           C          CA        25         1
#>    [4]     SSR3       176      - |           A          AG         8        16
#>    [5]     SSR3       198      - |           A          AG        24         1
#>    ...      ...       ...    ... .         ...         ...       ...       ...
#>   [70]     DHFR       397      - |           A          AT        31         1
#>   [71]     DHFR       399      - |           G          GA        28         1
#>   [72]     DHFR       423      - |           T          TC        31         1
#>   [73]     DHFR       430      - |           A           -        33         0
#>   [74]     DHFR       513      - |           T           -        21         0
#>               nA        nT        nC        nG        nN        nX
#>        <integer> <integer> <integer> <integer> <integer> <integer>
#>    [1]         0        15         0         0         0         7
#>    [2]         0         0        21         1         0         0
#>    [3]         1         0        25         0         0         0
#>    [4]         8         0         0        16         0         0
#>    [5]        24         0         0         1         0         1
#>    ...       ...       ...       ...       ...       ...       ...
#>   [70]        31         1         0         0         0         2
#>   [71]         1         0         0        28         0         6
#>   [72]         0        31         1         0         0         2
#>   [73]        33         0         0         0         0         2
#>   [74]         0        21         0         0         0         0
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
#> 
#> [[2]]
#> GRanges object with 74 ranges and 10 metadata columns:
#>        seqnames    ranges strand |         Ref         Var      nRef      nVar
#>           <Rle> <IRanges>  <Rle> | <character> <character> <integer> <integer>
#>    [1]     SSR3       102      - |           T          TG        12         1
#>    [2]     SSR3       125      - |           C           -        17         0
#>    [3]     SSR3       156      - |           C           -        16         0
#>    [4]     SSR3       176      - |           A           -        16         0
#>    [5]     SSR3       198      - |           A           -        15         0
#>    ...      ...       ...    ... .         ...         ...       ...       ...
#>   [70]     DHFR       397      - |           A           -        43         0
#>   [71]     DHFR       399      - |           G           -        43         0
#>   [72]     DHFR       423      - |           T           -        42         0
#>   [73]     DHFR       430      - |           A          AT        37         1
#>   [74]     DHFR       513      - |           T          TC        35         1
#>               nA        nT        nC        nG        nN        nX
#>        <integer> <integer> <integer> <integer> <integer> <integer>
#>    [1]         0        12         0         1         0         6
#>    [2]         0         0        17         0         0         1
#>    [3]         0         0        16         0         0         1
#>    [4]        16         0         0         0         0         0
#>    [5]        15         0         0         0         0         2
#>    ...       ...       ...       ...       ...       ...       ...
#>   [70]        43         0         0         0         0         2
#>   [71]         0         0         0        43         0         3
#>   [72]         0        42         0         0         0         1
#>   [73]        37         1         0         0         0         5
#>   [74]         0        35         1         0         0         0
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```

To facilitate comparisons across groups, the pileups can be stored in a
`RangedSummarizedExperiment`.

``` r
create_se(plps)
#> class: RangedSummarizedExperiment 
#> dim: 74 2 
#> metadata(0):
#> assays(7): Var nRef ... nC nG
#> rownames(74): SSR3_102_2 SSR3_125_2 ... DHFR_430_2 DHFR_513_2
#> rowData names(1): Ref
#> colnames(2): sample_1 sample_2
#> colData names(1): sample
```

## Related work

The functionality in `raer` builds off of previously published methods
and software:

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
