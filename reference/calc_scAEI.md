# Calculate the Adenosine Editing Index (AEI) in single cells

The Adenosine Editing Index describes the magnitude of A-to-I editing in
a sample. The index is a weighted average of editing events (G bases)
observed at A positions. The vast majority A-to-I editing occurs in ALU
elements in the human genome, and these regions have a high A-to-I
editing signal compared to other regions such as coding exons. This
function will examine enumerate edited and non-edited base counts at the
supplied sites and return summary AEI metric per cell. Potential editing
sites within repeat regions can be generated using `get_scAEI_sites()`.

## Usage

``` r
calc_scAEI(
  bamfiles,
  sites,
  cell_barcodes,
  param = FilterParam(),
  edit_from = "A",
  edit_to = "G",
  output_dir = NULL,
  return_sce = FALSE,
  ...
)

get_scAEI_sites(fasta, genes, alus, edit_from = "A", edit_to = "G")
```

## Arguments

- bamfiles:

  a path to a BAM file (for 10x libraries), or a vector of paths to BAM
  files (smart-seq2). Can be supplied as a character vector, `BamFile`,
  or `BamFileList`.

- sites:

  a GRanges object produced by `get_scAEI_sites()` containing sites to
  process.

- cell_barcodes:

  A character vector of single cell barcodes to process. If processing
  multiple BAM files (e.g. smart-seq-2), provide a character vector of
  unique identifiers for each input BAM, to name each BAM file in the
  output files.

- param:

  object of class
  [`FilterParam()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)
  which specify various filters to apply to reads and sites during
  pileup.

- edit_from:

  This should correspond to the base (`A`, `C`, `G`, `T`) you expect in
  the reference. Ex. for A to I editing events, this would be `A`.

- edit_to:

  This should correspond to the base (`A`, `C`, `G`, `T`) you expect in
  an edited site. Ex. for A to I editing events, this would be `G`.

- output_dir:

  Output directory for `nRef` and `nAlt` sparseMatrix files. If NULL, a
  temporary directory will be used.

- return_sce:

  if `TRUE`, data is returned as a SingleCellExperiment, if `FALSE` a
  `DataFrame` containing computed AEI values will be returned.

- ...:

  additional arguments to
  [`pileup_cells()`](https://rnabioco.github.io/raer/reference/pileup_cells.md)

- fasta:

  Path to a genome fasta file

- genes:

  A `GRanges` object with gene coordinates.Alternatively a `TxDb`
  object, which if supplied, will be used extract gene coordinates.

- alus:

  `GRanges` with repeat regions to query for calculating the AEI,
  typically ALU repeats. The strand of the supplied intervals will be
  ignored for defining repeat regions.

## Value

A `DataFrame` containing computed `AEI` values, count of editing events
(`n_alt`), and count of reference events (`n_ref`) per cell. If
`return_sce` is `TRUE`, then a `SingleCellExperiment` is returned with
the AEI values stored in the `colData`.

## References

Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide quantification of
ADAR adenosine-to-inosine RNA editing activity. Nat Methods 16,
1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9

## Examples

``` r
suppressPackageStartupMessages(library(Rsamtools))
library(GenomicRanges)

bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
bai <- indexBam(bam_fn)

# cell barcodes to query
cbs <- c("TGTTTGTTCCATCCGT-1", "CAACCAACATAATCGC-1", "TGGAACTCAAGCTGTT-1")

# genes used to infer transcribed strand
genes_gr <- GRanges(c(
    "2:100-400:-",
    "2:500-605:-",
    "2:600-680:+"
))

# alu intervals
alus_gr <- GRanges(c(
    "2:110-380",
    "2:510-600",
    "2:610-670"
))

# genome fasta file, used to find A bases
fa_fn <- raer_example("mouse_tiny.fasta")

# get positions of potential A -> G changes in alus
sites <- get_scAEI_sites(fa_fn, genes_gr, alus_gr)

fp <- FilterParam(
    library_type = "fr-second-strand",
    min_mapq = 255
)
calc_scAEI(bam_fn, sites, cbs, fp)
#> DataFrame with 3 rows and 3 columns
#>                          AEI     n_alt     n_ref
#>                    <numeric> <numeric> <numeric>
#> TGTTTGTTCCATCCGT-1   5.08475         3        56
#> CAACCAACATAATCGC-1  20.00000         3        12
#> TGGAACTCAAGCTGTT-1   0.00000         0         9
```
