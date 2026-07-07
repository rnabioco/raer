# Package index

## Pileup

- [`pileup_cells()`](https://rnabioco.github.io/raer/reference/pileup_cells.md)
  : Generate base counts per cell
- [`pileup_sites()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)
  [`FilterParam()`](https://rnabioco.github.io/raer/reference/pileup_sites.md)
  : Generate base counts using pileup

## Calculations

- [`calc_AEI()`](https://rnabioco.github.io/raer/reference/calc_AEI.md)
  : Calculate the Adenosine Editing Index (AEI)
- [`calc_confidence()`](https://rnabioco.github.io/raer/reference/calc_confidence.md)
  : Calculate confidence score for observing editing
- [`calc_edit_frequency()`](https://rnabioco.github.io/raer/reference/calc_edit_frequency.md)
  : Adds editing frequencies
- [`calc_scAEI()`](https://rnabioco.github.io/raer/reference/calc_scAEI.md)
  [`get_scAEI_sites()`](https://rnabioco.github.io/raer/reference/calc_scAEI.md)
  : Calculate the Adenosine Editing Index (AEI) in single cells

## Filtering

- [`filter_clustered_variants()`](https://rnabioco.github.io/raer/reference/filter_clustered_variants.md)
  : Filter out clustered sequence variants
- [`filter_multiallelic()`](https://rnabioco.github.io/raer/reference/filter_multiallelic.md)
  : Filter out multi-allelic sites
- [`filter_splice_variants()`](https://rnabioco.github.io/raer/reference/filter_splice_variants.md)
  : Filter out sites near splice sites

## Annotation

- [`annot_from_gr()`](https://rnabioco.github.io/raer/reference/annot_from_gr.md)
  : Annotate sites using GRanges object
- [`annot_snps()`](https://rnabioco.github.io/raer/reference/annot_snps.md)
  : Annotate known SNP positions
- [`get_overlapping_snps()`](https://rnabioco.github.io/raer/reference/get_overlapping_snps.md)
  : Retrieve SNPs overlapping intervals
- [`get_splice_sites()`](https://rnabioco.github.io/raer/reference/get_splice_sites.md)
  : Extract regions surrounding splice sites

## Bulk editing analysis

- [`find_de_sites()`](https://rnabioco.github.io/raer/reference/find_de_sites.md)
  : Perform differential editing
- [`make_de_object()`](https://rnabioco.github.io/raer/reference/make_de_object.md)
  : Make summarized experiment object for differential editing analysis

## Single cell editing analysis

- [`find_scde_sites()`](https://rnabioco.github.io/raer/reference/find_scde_sites.md)
  : Identify sites with differential editing between cells in single
  cell datasets
- [`find_mispriming_sites()`](https://rnabioco.github.io/raer/reference/find_mispriming_sites.md)
  : Find regions with oligodT mispriming

## Utilities

- [`correct_strand()`](https://rnabioco.github.io/raer/reference/correct_strand.md)
  : Apply strand correction using gene annotations
- [`calc_scAEI()`](https://rnabioco.github.io/raer/reference/calc_scAEI.md)
  [`get_scAEI_sites()`](https://rnabioco.github.io/raer/reference/calc_scAEI.md)
  : Calculate the Adenosine Editing Index (AEI) in single cells
- [`mock_rse()`](https://rnabioco.github.io/raer/reference/mock_rse.md)
  : Generate a small RangedSummarizedExperiment object for tests and
  examples
- [`raer_example()`](https://rnabioco.github.io/raer/reference/raer_example.md)
  : Provide working directory for raer example files.
- [`read_sparray()`](https://rnabioco.github.io/raer/reference/read_sparray.md)
  : Read sparseMatrix produced by pileup_cells()
