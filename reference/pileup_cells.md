# Generate base counts per cell

This function processes scRNA-seq library to enumerate base counts for
Reference (unedited) or Alternate (edited) bases at specified sites in
single cells. `pileup_cells` can process droplet scRNA-seq libraries,
from a BAM file containing a cell-barcode and UMI, or well-based
libraries that do not contain cell-barcodes.

The `sites` parameter specifies sites to quantify. This must be a
[GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with 1 base intervals, a strand (+ or -), and supplemented with
metadata columns named `REF` and `ALT` containing the reference and
alternate base to query. See examples for the required format.

At each site, bases from overlapping reads will be examined, and counts
of each ref and alt base enumerated for each cell-barcode present. A
single base will be counted once for each UMI sequence present in each
cell.

## Usage

``` r
pileup_cells(
  bamfiles,
  sites,
  cell_barcodes,
  output_directory,
  chroms = NULL,
  umi_tag = "UB",
  cb_tag = "CB",
  param = FilterParam(),
  BPPARAM = SerialParam(),
  return_sce = TRUE,
  verbose = FALSE
)
```

## Arguments

- bamfiles:

  a path to a BAM file (for droplet scRNA-seq), or a vector of paths to
  BAM files (Smart-seq2). Can be supplied as a character vector,
  [BamFile](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html), or
  [BamFileList](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html).

- sites:

  a GRanges object containing sites to process. See examples for valid
  formatting.

- cell_barcodes:

  A character vector of single cell barcodes to process. If processing
  multiple BAM files (e.g. Smart-seq2), provide a character vector of
  unique identifiers for each input BAM, to name each BAM file in the
  output files.

- output_directory:

  Output directory for output matrix files. The directory will be
  generated if it doesn't exist.

- chroms:

  A character vector of chromosomes to process. If supplied, only sites
  present in the listed chromosomes will be processed

- umi_tag:

  tag in BAM containing the UMI sequence

- cb_tag:

  tag in BAM containing the cell-barcode sequence

- param:

  object of class
  [`FilterParam()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)
  which specify various filters to apply to reads and sites during
  pileup. Note that the `min_depth` and `min_variant_reads` parameters
  if set \> 0 specify the number of reads from any cell required in
  order to report a site. E.g. if `min_variant_reads` is set to 2, then
  at least 2 reads (from any cell) must have a variant in order to
  report the site. Setting `min_depth` and `min_variant_reads` to 0
  reports all sites present in the `sites` object. The following options
  are not enabled for pileup_cells(): `max_mismatch_type`,
  `homopolymer_len`, and `min_allelic_freq`.

- BPPARAM:

  BiocParallel instance. Parallel computation occurs across chromosomes.

- return_sce:

  if `TRUE`, data is returned as a SingleCellExperiment, if `FALSE` a
  character vector of the output files, specified by `outfile_prefix`,
  will be returned.

- verbose:

  Display messages

## Value

Returns either a SingleCellExperiment or character vector of paths to
the sparseMatrix files produced. The SingleCellExperiment object is
populated with two assays, `nRef` and `nAlt`, which represent base
counts for the reference and alternate alleles. The `rowRanges()` will
contain the genomic interval for each site, along with `REF` and `ALT`
columns. The rownames will be populated with the format
`site_[seqnames]_[position(1-based)]_[strand]_[allele]`, with `strand`
being encoded as 1 = +, 2 = -, and 3 = \*, and allele being `REF` +
`ALT`.

If `return_sce` is `FALSE` then a character vector of paths to the
sparseMatrix files (`barcodes.txt.gz`, `sites.txt.gz`, `counts.mtx.gz`),
will be returned. These files can be imported using
[`read_sparray()`](https://rnabioco.github.io/raer/reference/read_sparray.md).

## See also

Other pileup:
[`pileup_sites()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)

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
sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
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

# example of processing multiple Smart-seq2 style libraries

many_small_bams <- rep(bam_fn, 10)
bam_ids <- LETTERS[1:10]

# for unstranded libraries, sites and alleles should be provided on + strand
gr <- GRanges(c("2:579:+", "2:625:+", "2:645:+", "2:589:+", "2:601:+"))
gr$REF <- c(rep("T", 4), "A")
gr$ALT <- c(rep("C", 4), "G")

fp <- FilterParam(
    library_type = "unstranded",
    remove_overlaps = TRUE
)

sce <- pileup_cells(many_small_bams,
    sites = gr,
    cell_barcodes = bam_ids,
    cb_tag = NULL,
    umi_tag = NULL,
    outdir,
    param = fp
)
sce
#> class: SingleCellExperiment 
#> dim: 4 10 
#> metadata(0):
#> assays(2): nRef nAlt
#> rownames(4): site_2_579_1_TC site_2_625_1_TC site_2_589_1_TC
#>   site_2_601_1_AG
#> rowData names(2): REF ALT
#> colnames(10): A B ... I J
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

unlink(bai)
```
