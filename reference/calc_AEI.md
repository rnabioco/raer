# Calculate the Adenosine Editing Index (AEI)

The Adenosine Editing Index describes the magnitude of A-to-I editing in
a sample. The index is a weighted average of editing events (G bases)
observed at A positions. The vast majority A-to-I editing occurs in ALU
elements in the human genome, and these regions have a high A-to-I
editing signal compared to other regions such as coding exons. This
function will perform pileup at specified repeat regions and return a
summary AEI metric.

## Usage

``` r
calc_AEI(
  bamfiles,
  fasta,
  alu_ranges = NULL,
  txdb = NULL,
  snp_db = NULL,
  param = FilterParam(),
  BPPARAM = SerialParam(),
  verbose = FALSE
)
```

## Arguments

- bamfiles:

  character vector of paths to indexed bam files. If a named character
  vector is supplied the names will be used in the output.

- fasta:

  fasta filename

- alu_ranges:

  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  with regions to query for calculating the AEI, typically ALU repeats.

- txdb:

  A TxDb object, if supplied, will be used to subset the alu_ranges to
  those found overlapping genes. Alternatively a
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object with gene coordinates. If the `library_type`, specified by
  `FilterParam`, is `unstranded` then the TxDb will be used to correct
  the strandness relative to the reference and is a required parameter.

- snp_db:

  either a
  [SNPlocs](https://rdrr.io/pkg/BSgenome/man/SNPlocs-class.html),
  [GPos](https://rdrr.io/pkg/GenomicRanges/man/GPos-class.html), or
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object. If supplied, will be used to exclude polymorphic positions
  prior to calculating the AEI. If `calc_AEI()` will be used many times,
  one will save time by first identifying SNPs that overlap the supplied
  `alu_ranges`, and passing these as a
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) to
  `snp_db` rather than supplying all known SNPs (see
  [`get_overlapping_snps()`](https://rnabioco.github.io/raer/reference/get_overlapping_snps.md)).

- param:

  object of class
  [`FilterParam()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)
  which specify various filters to apply to reads and sites during
  pileup.

- BPPARAM:

  A BiocParallelParam object for specifying parallel options for
  operating over chromosomes.

- verbose:

  report progress on each chromosome?

## Value

A named list containing:

- `AEI`: a matrix of AEI index values computed for all allelic
  combinations, one row for each supplied bam file.

- `AEI_per_chrom`: a data.frame containing values computed for each
  chromosome

## References

Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide quantification of
ADAR adenosine-to-inosine RNA editing activity. Nat Methods 16,
1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9

## Examples

``` r
suppressPackageStartupMessages(library(Rsamtools))

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
bams <- c(bamfn, bam2fn)
names(bams) <- c("ADAR1KO", "WT")

fafn <- raer_example("human.fasta")
mock_alu_ranges <- scanFaIndex(fafn)

res <- calc_AEI(bams, fafn, mock_alu_ranges)
res$AEI
#>                A_C       A_G       A_T       C_A        C_G       C_T       G_A
#> ADAR1KO 0.00927902 0.1019746 0.0278319 0.1237964 0.00000000 0.1100564 0.2364360
#> WT      0.00000000 2.2985864 0.1292293 0.2142033 0.01650982 0.2635046 0.3170769
#>                G_C        G_T       T_A        T_C        T_G
#> ADAR1KO 0.01247194 0.04986909 0.0108272 0.03247456 0.03247456
#> WT      0.00000000 0.00000000 0.0000000 0.01161980 0.01161980
```
