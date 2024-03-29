
# raer 1.1.3

* Fixed bug in smart-seq2 pileup_cells reporting

# raer 1.0.0

* raer added to the Bioconductor project

# raer 0.99.12 

* Changes made to prepare for bioc submission

# raer 0.99.11

* Replaced base R fisher test with c-wrapper to call htslib fisher test, which speeds up execution many fold. 

# raer 0.99.10

* The options to write to tabix indexed output files have been removed from `pileup_sites()` as they have limited utility and introduce unwanted code complexity. 

# raer 0.99.9

* The `genomic-unstranded` option for the `library-type` argument in `FilterParam()` has been renamed to `unstranded`, and the `unstranded` option has been removed. 

# raer 0.99.8

* Function arguments involving a fasta file have been renamed to all be `fasta`

# raer 0.99.7

* added a single cell specific AEI calculation (`calc_scAEI()`)

# raer 0.99.6 

* added method to count base consensus base when counting UMIs with `pileup_cells()` using the sum of base qualities to select consensus. 

# raer 0.99.5

* `pileup_cells()` now allows for multiple alleles to be queried at a site.

* Fixed an indexing bug in `pileup_cells()` that misassigned sites to counts.

# raer 0.99.4

* `annot_snps` will now compare the site allele to the SNP allele and return a 
column `snp_matches_site` indicating if the site matches the snp. 

* added new function, `find_scde_sites()` to identify differentially editing 
sites in single cell data using fishers exact tests. 

* `pileup_cells` now respects the `min_depth` and `min_variant_reads` FilterParameters.   

# raer 0.99.3

* support `BamFile` and `BamFileList` inputs to `pileup_sites()` and `pileup_cells()`, which provides an option to provide custom BAI index file names.

# raer 0.99.2

* rename `prep_for_de()` and `perform_de()` to `make_de_object()` and `find_de_sites()`. 

# raer 0.99.1 

* default values for `edit_from` and `edit_to` for `calc_edit_frequency()` have
been changed to `A` and `G` respectively.

* renamed `type` argument in `perform_de` to `test` and removed `type` argument in
`prep_for_de`

# raer 0.99.0

* added support for processing multiple BAM files with `calc_AEI()`. 

* Dropped minimally used `bad_reads` and `reads` parameters from `pileup_sites()`

* Added utility to screen scRNA-seq bam files for regions with oligo-dT mispriming (`find_mispriming_sites()`).

* add option to query ref and alt SNP alleles 

* added tests for SummarizedExperiment filtering approaches 

* added a strand bias stat `sor` using approach from GATK (StrandOddsRatio), and
  a confidence score `calc_confidence()` from SAILOR pipeline.

* 'N' bases in read or reference are ignored 

* Removed outdated or unused functionality:
  - bed indexing (`indexBed` and related C code)
  - bam tag indexing (`build_tag_index`, `show_tag_index`, `get_tag_bam`, )
  - bam tag index based single cell approach (`sc_editing`)
  - bam tag indexing C code from bri (`src/bri/*`)
  - sparse matrix merging for `merge_pileups()`.
  - unneeded utilities (`filter_by_coverage`)
  - Remaining (and mostly unused) Rcpp code
  - Removed fastmap, Rcpp, zlibbioc, RColorBrewer, and BiocGenerics dependencies
  - Removed system requirements for C libraries used by bri
  
* The bed indexing used in `pileup_sites()` has been replaced with the region indexing approach from `pileup_cells()`. 

* `pileup_sites()` now requires a GRanges object rather than a bed file. The `bedfile `parameter has been removed and replaced with a `sites` parameter.  

* Renamed `Ref` and `Var` output columns to `REF` and `ALT` and `nVar` was renamed to `nAlt`. This provides consistency with VCF format and consistency across `pileup_cells()` and `pileup_sites()` function calls

* `pileup_cells()` gained functionality to process multiple smart-seq2 style bam files.

* Changed `filterParam` argument in `pileup_sites` and `pileup_cells` to `param` for simplicity.

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
