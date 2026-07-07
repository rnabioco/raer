# Identify sites with differential editing between cells in single cell datasets

Compare editing frequencies between clusters or celltypes. REF and ALT
counts from each cluster are pooled to create pseudobulk estimates. Each
pair of clusters are compared using fisher exact tests. Statistics are
aggregated across each pairwise comparison using
[scran::combineMarkers](https://rdrr.io/pkg/scran/man/combineMarkers.html).

## Usage

``` r
find_scde_sites(sce, group, rowData = FALSE, BPPARAM = SerialParam(), ...)
```

## Arguments

- sce:

  SingleCellExperiment object with `nRef` and `nAlt` assays.

- group:

  column name from colData used to define groups to compare.

- rowData:

  if TRUE,
  [rowData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  from the input SingleCellExperiment will be included in the output
  DataFrames

- BPPARAM:

  BiocParallel backend for control how parallel computations are
  performed.

- ...:

  Additional arguments passed to
  [scran::combineMarkers](https://rdrr.io/pkg/scran/man/combineMarkers.html)

## Value

A named list of
[DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)s
containing results for each cluster specified by `group`. The difference
in editing frequencies between cluster pairs are denoted as `dEF`. See
[scran::combineMarkers](https://rdrr.io/pkg/scran/man/combineMarkers.html)
for a description of additional output fields.

## Examples

``` r

### generate example data ###

library(Rsamtools)
library(GenomicRanges)
bam_fn <- raer_example("5k_neuron_mouse_possort.bam")

gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
gr$REF <- c(rep("A", 4), "T")
gr$ALT <- c(rep("G", 4), "C")

cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
cbs <- na.omit(cbs)

outdir <- tempdir()
bai <- indexBam(bam_fn)

fp <- FilterParam(library_type = "fr-second-strand")
sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)

# mock some clusters
set.seed(42)
sce$clusters <- paste0("cluster_", sample(1:3, ncol(sce), replace = TRUE))
res <- find_scde_sites(sce, "clusters")
#> → 244 cells had no REF or ALT counts and were excluded from the analysis
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'summarizeAssayByGroup' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' or 'beachmat::tatami.sums.by.group' instead.
#> Warning: 'summarizeAssayByGroup' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' or 'beachmat::tatami.sums.by.group' instead.
#> Warning: 'scran::combineMarkers' is deprecated.
#> Use 'scrapper::summarizeEffects' instead.
#> See help("Deprecated")
#> Warning: 'summaryMarkerStats' is deprecated.
#> See help("Deprecated")
#> Warning: 'sumCountsAcrossCells' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' instead.
#> See help("Deprecated")
#> Warning: 'summarizeAssayByGroup' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' or 'beachmat::tatami.sums.by.group' instead.
#> Warning: 'numDetectedAcrossCells' is deprecated.
#> Use 'scrapper::computeRnaQcMetrics' instead.
#> See help("Deprecated")
#> Warning: 'summarizeAssayByGroup' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' or 'beachmat::tatami.sums.by.group' instead.
res[[1]]
#> DataFrame with 4 rows and 10 columns
#>                 self.average other.average self.detected other.detected
#>                    <numeric>     <numeric>     <numeric>      <numeric>
#> site_2_579_2_AG     0.690208      0.653495      0.727273       0.707399
#> site_2_625_2_AG     0.926220      0.991736      0.838384       0.901309
#> site_2_589_2_AG     0.862085      0.831895      0.868687       0.861512
#> site_2_601_2_TC     0.838316      0.850864      0.858586       0.871928
#>                       Top   p.value       FDR summary.dEF dEF.cluster_2
#>                 <integer> <numeric> <numeric>   <numeric>     <numeric>
#> site_2_579_2_AG         1  0.872535         1  0.01784427    0.02257960
#> site_2_625_2_AG         1  1.000000         1 -0.03537552   -0.03537552
#> site_2_589_2_AG         3  1.000000         1 -0.00772495   -0.00772495
#> site_2_601_2_TC         4  1.000000         1 -0.00121012   -0.00121012
#>                 dEF.cluster_3
#>                     <numeric>
#> site_2_579_2_AG   0.017844268
#> site_2_625_2_AG  -0.008314991
#> site_2_589_2_AG  -0.000749251
#> site_2_601_2_TC   0.000124844
```
