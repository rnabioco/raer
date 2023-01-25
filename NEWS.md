# raer 0.99.0

* Rename `remove_clustered_variants`, `remove_multiallelic`, and `remove_splice_variants` 
  to `filter_*` for consistency.

* Rewrote and renamed the single cell editing function `sc_editing` to `pileup_cells()`. `pileup_cells()` does not require sorting and index by cell barcode, uses a new format to specify sites to query and requires providing the reference and alternate alleles of interest, writes to disk in a sparse matrix compatible format to reduce memory usage, and should have more performance as there is no need to query a fasta index. 

* Implemented method to collapse reads with duplicate UMIs.

* Added option to filter sites in pileup based on number of reads containing a variant (#54)

* Added a `NEWS.md` file to track changes to the package.
