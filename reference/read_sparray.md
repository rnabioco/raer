# Read sparseMatrix produced by pileup_cells()

Read in tables produced by
[`pileup_cells()`](https://rnabioco.github.io/raer/reference/pileup_cells.md)
which are an extension of the matrixMarket sparse matrix format to store
values for more than 1 matrix.

The .mtx.gz files are formatted with columns:

1.  row index (0 based)

2.  column index (0 based)

3.  values for sparseMatrix \#1 (nRef)

4.  values for sparseMatrix \#2 (nAlt)

## Usage

``` r
read_sparray(mtx_fn, sites_fn, bc_fn, site_format = c("coordinate", "index"))
```

## Arguments

- mtx_fn:

  .mtx.gz file path

- sites_fn:

  sites.txt.gz file path

- bc_fn:

  bcs.txt.gz file path

- site_format:

  one of `coordinate` or `index`, `coordinate` will populate a
  SingleCellExperiment with rowRanges and rownames corresponing to
  genomic intervals, whereas \`index“ will only add row indices to the
  rownames.

## Value

a `SingleCellExperiment` object populated with `nRef` and `nAlt` assays.

## Examples

``` r
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
mtx_fns <- pileup_cells(bam_fn, gr, cbs, outdir, return_sce = FALSE)
sce <- read_sparray(mtx_fns[1], mtx_fns[2], mtx_fns[3])
sce
#> class: SingleCellExperiment 
#> dim: 4 556 
#> metadata(0):
#> assays(2): nRef nAlt
#> rownames(4): site_2_579_2_AG site_2_625_2_AG site_2_589_2_AG
#>   site_2_601_2_TC
#> rowData names(2): REF ALT
#> colnames(556): TGGAACTCAAGCTGTT-1 TACTTCAGTAACCCTA-1 ...
#>   TGTACAGTCTTCGTGC-1 TGTTGAGGTGACTGAG-1
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

unlink(bai)
```
