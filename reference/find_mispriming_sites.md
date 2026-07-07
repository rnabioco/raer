# Find regions with oligodT mispriming

OligodT will prime at A-rich regions in an RNA. Reverse transcription
from these internal priming sites will install an oligodT sequence at
the 3' end of the cDNA. Sequence variants within these internal priming
sites are enriched for variants converting the genomic sequence to the A
encoded by the oligodT primer. Trimming poly(A) from the 3' ends of
reads reduces but does not eliminate these signals

This function will identify regions that are enriched for mispriming
events. Reads that were trimmed to remove poly(A) (encoded in the pa tag
by 10x Genomics) are identified. The aligned 3' positions of these reads
are counted, and sites passing thresholds (at least 2 reads) are
retained as possible sites of mispriming. Be default regions 5 bases
upstream and 20 bases downstream of these putative mispriming sites are
returned.

## Usage

``` r
find_mispriming_sites(
  bamfile,
  fasta,
  pos_5p = 5,
  pos_3p = 20,
  min_reads = 2,
  tag = "pa",
  tag_values = 3:300,
  n_reads_per_chunk = 1e+06,
  verbose = TRUE
)
```

## Arguments

- bamfile:

  path to bamfile

- fasta:

  path to fasta file

- pos_5p:

  distance 5' of mispriming site to define mispriming region

- pos_3p:

  distance 3' of mispriming site to define mispriming region

- min_reads:

  minimum required number of reads at a mispriming site

- tag:

  bam tag containing number of poly(A) bases trimmed

- tag_values:

  range of values required for read to be considered

- n_reads_per_chunk:

  number of reads to process in memory, see
  [`Rsamtools::BamFile()`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html)

- verbose:

  if true report progress

## Value

A GenomicsRanges containing regions enriched for putative mispriming
events. The `n_reads` column specifies the number of polyA trimmed reads
overlapping the mispriming region. `mean_pal` indicates the mean length
of polyA sequence trimmed from reads overlapping the region. The
`n_regions` column specifies the number overlapping independent regions
found in each chunk (dictated by `n_reads_per_chunk`). The `A_freq`
column indicates the frequency of A bases within the region.

## Examples

``` r
bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
fa_fn <- raer_example("mouse_tiny.fasta")
find_mispriming_sites(bam_fn, fa_fn)
#> working on 2:580-633:- to 2:580-633:-
#> GRanges object with 1 range and 4 metadata columns:
#>       seqnames    ranges strand |  mean_pal   n_reads n_regions    A_freq
#>          <Rle> <IRanges>  <Rle> | <numeric> <integer> <integer> <numeric>
#>   [1]        2   560-585      - |   34.6667         6         1  0.615385
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome
```
