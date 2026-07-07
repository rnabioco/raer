# Annotate sites using GRanges object

Utility function to map annotations from GRanges to rowData of
SummarizedExperiment or to mcols of GRanges object. If multiple features
overlap then they will be concatenated with the specified separtor
string.

## Usage

``` r
annot_from_gr(obj, gr, cols_to_map, RLE = TRUE, sep = ",", ...)
```

## Arguments

- obj:

  RangedSummarizedExperiment or GRanges object

- gr:

  GRanges with annotations to map to obj

- cols_to_map:

  character vector of columns from GRanges to map to
  SummarizedExperiment. If the vector has names, the names will be the
  column names in the output.

- RLE:

  If TRUE, columns added will returned as
  [`S4Vectors::Rle()`](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html)
  vectors to reduce memory

- sep:

  separator string, defaults to comma.

- ...:

  additional arguments to pass to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)

## Value

Either a SummarizedExperiment or GRanges object with additional
annotations provided by the supplied GRanges object.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
rse_adar_ifn <- mock_rse()
gr <- GRanges(rep(c("SSR3", "SPCS3"), c(5, 15)),
    IRanges(seq(1, 500, by = 25), width = 50),
    strand = "+"
)

gr$feature <- sample(1:100, size = 20)
gr$id <- sample(LETTERS, size = 20)

rse <- annot_from_gr(rse_adar_ifn, gr, c(feature_set = "feature", "id"))
rowData(rse)
#> DataFrame with 74 rows and 6 columns
#>                         REF       rpbz        vdb       sor feature_set    id
#>                 <character>  <numeric>  <numeric> <numeric>       <Rle> <Rle>
#> site_SSR3_102_2           T   1.489645        Inf  1.403164          NA    NA
#> site_SSR3_125_2           C   0.356711        Inf  0.165499          NA    NA
#> site_SSR3_156_2           C   1.073919        Inf  1.442924          NA    NA
#> site_SSR3_176_2           A  -0.387238 0.00686469  1.987570          NA    NA
#> site_SSR3_198_2           A   1.040581        Inf  1.483784          NA    NA
#> ...                     ...        ...        ...       ...         ...   ...
#> site_DHFR_397_2           A -1.5715051        Inf   1.39896          NA    NA
#> site_DHFR_399_2           G -0.1203878        Inf   0.09602          NA    NA
#> site_DHFR_423_2           T -0.0468703        Inf   1.38985          NA    NA
#> site_DHFR_430_2           A -1.5389404        Inf   1.39019          NA    NA
#> site_DHFR_513_2           T -0.7160074        Inf   1.38637          NA    NA
```
