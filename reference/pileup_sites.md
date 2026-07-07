# Generate base counts using pileup

This function uses a pileup routine to examine numerate base counts from
alignments at specified sites, regions, or across all read alignments,
from one or more BAM files. Alignment and site filtering options are
controlled by the `FilterParam` class. A
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
object is returned, populated with base count statistics for each
supplied BAM file.

## Usage

``` r
pileup_sites(
  bamfiles,
  fasta,
  sites = NULL,
  region = NULL,
  chroms = NULL,
  param = FilterParam(),
  BPPARAM = SerialParam(),
  umi_tag = NULL,
  verbose = FALSE
)

FilterParam(
  max_depth = 10000,
  min_depth = 1L,
  min_base_quality = 20L,
  min_mapq = 0L,
  library_type = "fr-first-strand",
  bam_flags = NULL,
  only_keep_variants = FALSE,
  trim_5p = 0L,
  trim_3p = 0L,
  ftrim_5p = 0,
  ftrim_3p = 0,
  indel_dist = 0L,
  splice_dist = 0L,
  min_splice_overhang = 0L,
  homopolymer_len = 0L,
  max_mismatch_type = c(0L, 0L),
  read_bqual = c(0, 0),
  min_variant_reads = 0L,
  min_allelic_freq = 0,
  report_multiallelic = TRUE,
  remove_overlaps = TRUE
)
```

## Arguments

- bamfiles:

  a character vector,
  [BamFile](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html) or
  [BamFileList](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html)
  indicating 1 or more BAM files to process. If named, the names will be
  included in the
  [colData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  of the
  [RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  as a `sample` column, otherwise the names will be taken from the
  basename of the BAM file.

- fasta:

  path to genome fasta file used for read alignment. Can be provided in
  compressed gzip or bgzip format.

- sites:

  a [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object containing regions or sites to process.

- region:

  samtools region query string (i.e. `chr1:100-1000`). Can be combined
  with sites, in which case sites will be filtered to keep only sites
  within the region.

- chroms:

  chromosomes to process, provided as a character vector. Not to be used
  with the region parameter.

- param:

  object of class `FilterParam()` which specify various filters to apply
  to reads and sites during pileup.

- BPPARAM:

  A BiocParallel class to control parallel execution. Parallel
  processing occurs per chromosome and is disabled when run on a single
  region.

- umi_tag:

  The BAM tag containing a UMI sequence. If supplied, multiple reads
  with the same UMI sequence will only be counted once per position.

- verbose:

  if TRUE, then report progress and warnings.

- max_depth:

  maximum read depth considered at each site

- min_depth:

  min read depth needed to report site

- min_base_quality:

  min base quality score to consider read for pileup

- min_mapq:

  minimum required MAPQ score. Values for each input BAM file can be
  provided as a vector.

- library_type:

  read orientation, one of `fr-first-strand`, `fr-second-strand`, and
  `unstranded`. Unstranded library type will be reported with variants
  w.r.t the + strand. Values for each input BAM file can be provided as
  a vector.

- bam_flags:

  bam flags to filter or keep, use
  [`Rsamtools::scanBamFlag()`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  to generate.

- only_keep_variants:

  if TRUE, then only variant sites will be reported (FALSE by default).
  Values for each input BAM file can be provided as a vector.

- trim_5p:

  Bases to trim from 5' end of read alignments

- trim_3p:

  Bases to trim from 3' end of read alignments

- ftrim_5p:

  Fraction of bases to trim from 5' end of read alignments

- ftrim_3p:

  Fraction of bases to trim from 3' end of read alignments

- indel_dist:

  Exclude read if site occurs within given distance from indel event in
  the read

- splice_dist:

  Exclude read if site occurs within given distance from splicing event
  in the read

- min_splice_overhang:

  Exclude read if site is located adjacent to splice site with an
  overhang less than given length.

- homopolymer_len:

  Exclude site if occurs within homopolymer of given length

- max_mismatch_type:

  Exclude read if it has X different mismatch types (e.g A-to-G, G-to-C,
  C-to-G, is 3 mismatch types) or Y \# of mismatches, must be supplied
  as a integer vector of length 2. e.g. c(X, Y).

- read_bqual:

  Exclude read if more than X percent of the bases have base qualities
  less than Y. Numeric vector of length 2. e.g. c(0.25, 20)

- min_variant_reads:

  Required number of reads containing a variant for a site to be
  reported. Calculated per bam file, such that if 1 bam file has \>=
  min_variant_reads, then the site will be reported.

- min_allelic_freq:

  minimum allelic frequency required for a variant to be reported in ALT
  assay.

- report_multiallelic:

  if TRUE, report sites with multiple variants passing filters. If
  FALSE, site will not be reported.

- remove_overlaps:

  if TRUE, enable read pair overlap detection, which will count only 1
  read in regions where read pairs overlap using the htslib algorithm.
  In brief for each overlapping base pair the base quality of the base
  with the lower quality is set to 0, which discards it from being
  counted.

## Value

A
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
object populated with multiple assays:

- `ALT`: Alternate base(s) found at each position

- `nRef`: \# of reads supporting the reference base

- `nAlt`: \# of reads supporting an alternate base

- `nA`: \# of reads with A

- `nT`: \# of reads with T

- `nC`: \# of reads with C

- `nG`: \# of reads with G

The `rowRanges()` contains the genomic interval for each site, along
with:

- `REF`: The reference base

- `rpbz`: Mann-Whitney U test of Read Position Bias from bcftools,
  extreme negative or positive values indicate more bias.

- `vdb`: Variant Distance Bias for filtering splice-site artifacts from
  bcftools, lower values indicate more bias.

- `sor` Strand Odds Ratio Score, strand bias estimated by the Symmetric
  Odds Ratio test, based on GATK code. Higher values indicate more bias.

The rownames will be populated with the format
`site_[seqnames]_[position(1-based)]_[strand]`, with `strand` being
encoded as 1 = +, 2 = -, and 3 = \*.

## See also

Other pileup:
[`pileup_cells()`](https://rnabioco.github.io/raer/reference/pileup_cells.md)

## Examples

``` r
library(SummarizedExperiment)
bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")

rse <- pileup_sites(bamfn, fafn)

fp <- FilterParam(only_keep_variants = TRUE, min_depth = 55)
pileup_sites(bamfn, fafn, param = fp)
#> class: RangedSummarizedExperiment 
#> dim: 7 1 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(7): site_DHFR_207_2 site_DHFR_208_2 ... site_DHFR_300_2
#>   site_DHFR_332_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(1): SRR5564269_Aligned.sortedByCoord.out.md.bam
#> colData names(1): sample


# using multiple bam files

bams <- rep(c(bamfn, bam2fn), each = 3)
sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
names(bams) <- sample_ids

fp <- FilterParam(only_keep_variants = TRUE)
rse <- pileup_sites(bams, fafn, param = fp)
rse
#> class: RangedSummarizedExperiment 
#> dim: 74 6 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(6): KO1 KO2 ... WT2 WT3
#> colData names(1): sample

rse$condition <- substr(rse$sample, 1, 2)
assays(rse)
#> List of length 7
#> names(7): ALT nRef nAlt nA nT nC nG

colData(rse)
#> DataFrame with 6 rows and 2 columns
#>          sample   condition
#>     <character> <character>
#> KO1         KO1          KO
#> KO2         KO2          KO
#> KO3         KO3          KO
#> WT1         WT1          WT
#> WT2         WT2          WT
#> WT3         WT3          WT

rowRanges(rse)
#> GRanges object with 74 ranges and 4 metadata columns:
#>                   seqnames    ranges strand |         REF       rpbz
#>                      <Rle> <IRanges>  <Rle> | <character>  <numeric>
#>   site_SSR3_102_2     SSR3       102      - |           T   2.611800
#>   site_SSR3_125_2     SSR3       125      - |           C   0.623238
#>   site_SSR3_156_2     SSR3       156      - |           C   1.875144
#>   site_SSR3_176_2     SSR3       176      - |           A  -0.676423
#>   site_SSR3_198_2     SSR3       198      - |           A   1.817678
#>               ...      ...       ...    ... .         ...        ...
#>   site_DHFR_397_2     DHFR       397      - |           A -2.7341600
#>   site_DHFR_399_2     DHFR       399      - |           G -0.2094944
#>   site_DHFR_423_2     DHFR       423      - |           T -0.0815516
#>   site_DHFR_430_2     DHFR       430      - |           A -2.6781858
#>   site_DHFR_513_2     DHFR       513      - |           T -1.2475213
#>                           vdb       sor
#>                     <numeric> <numeric>
#>   site_SSR3_102_2 2.21621e-02  2.776147
#>   site_SSR3_125_2 2.21621e-02  0.478641
#>   site_SSR3_156_2 2.21621e-02  2.786335
#>   site_SSR3_176_2 2.00523e-05  2.154124
#>   site_SSR3_198_2 2.21621e-02  2.797103
#>               ...         ...       ...
#>   site_DHFR_397_2   0.0221621  2.775538
#>   site_DHFR_399_2   0.0221621  0.306541
#>   site_DHFR_423_2   0.0221621  2.773352
#>   site_DHFR_430_2   0.0221621  2.773426
#>   site_DHFR_513_2   0.0221621  2.772591
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome

# specifying regions to query using GRanges object

sites <- rowRanges(rse)
rse <- pileup_sites(bams, fafn, sites = sites)
rse
#> class: RangedSummarizedExperiment 
#> dim: 74 6 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(74): site_SSR3_102_2 site_SSR3_125_2 ... site_DHFR_430_2
#>   site_DHFR_513_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(6): KO1 KO2 ... WT2 WT3
#> colData names(1): sample

rse <- pileup_sites(bams, fafn, chroms = c("SPCS3", "DHFR"))
rse
#> class: RangedSummarizedExperiment 
#> dim: 1166 6 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(1166): site_SPCS3_1_1 site_SPCS3_2_1 ... site_DHFR_517_2
#>   site_DHFR_518_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(6): KO1 KO2 ... WT2 WT3
#> colData names(1): sample

rse <- pileup_sites(bams, fafn, region = "DHFR:100-101")
rse
#> class: RangedSummarizedExperiment 
#> dim: 2 6 
#> metadata(0):
#> assays(7): ALT nRef ... nC nG
#> rownames(2): site_DHFR_100_2 site_DHFR_101_2
#> rowData names(4): REF rpbz vdb sor
#> colnames(6): KO1 KO2 ... WT2 WT3
#> colData names(1): sample
```
