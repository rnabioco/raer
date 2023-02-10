# raer 0.99.0

* `pileup_cells()` gained functionality to process multiple smart-seq2 style bam files.

* Changed `filterParam` argument in `pileup_sites` and `pileup_cells` to `param` for
simplicity

* Added `FilterParam` to exclude multi-allelic sites `report_multiallelic`, or exclude reporting a variant in the Var assay based on allelic frequency (`min_allelic_freq`).

* The `bam_flags` parameter used in `pileup_sites` and `pileup_cells` has been moved into the `FilterParam` class. 

* The `bedindex` parameter for `pileup_sites` has been removed. This option is not needed
at the user level and is planned to be replaced by the regional indexing used in `pileup_cells()`.

* Added `FilterParam` option to trim reads based on fractional distance from 5' (`ftrim_5p`) or 3' end (`ftrim_3p`).

* Incorporated RBPZ and VDB statistics from bcftools, now returned as rowData columns 
when calling `pileup_sites`.

* A `RangedSummarizedExperiment` object is now directly returned from `pileup_sites`. Using `merge_pileups` is no longer necessary and is not an exported function. 

* Renamed `get_pileup` to `pileup_sites` and `create_se` to `merge_pileups`

* Rename `remove_clustered_variants`, `remove_multiallelic`, and `remove_splice_variants` 
  to `filter_*` for consistency.

* Rewrote and renamed the single cell editing function `sc_editing` to `pileup_cells()`. `pileup_cells()` does not require sorting and index by cell barcode, uses a new format to specify sites to query and requires providing the reference and alternate alleles of interest, writes to disk in a sparse matrix compatible format to reduce memory usage, and should have more performance as there is no need to query a fasta index. 

* Implemented method to collapse reads with duplicate UMIs.

* Added option to filter sites in pileup based on number of reads containing a variant (#54)

* Added a `NEWS.md` file to track changes to the package.
