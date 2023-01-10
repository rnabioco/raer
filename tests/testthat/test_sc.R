library(BiocParallel)
library(rtracklayer)
library(GenomicRanges)
library(SummarizedExperiment)

bamfn <- raer_example("5k_neuron_mouse_xf25_1pct_cbsort.bam")
idxfn <- build_tag_index(bamfn)
cbs <- show_tag_index(bamfn)$tag

cb_lst <- split(cbs, cut(seq_along(cbs), breaks = 5))
names(cb_lst) <- paste0("cluster", 1:5)

test_that("basic functionality works", {
  fp <- FilterParam(library_type = "fr-second-strand")
  se <- sc_editing(
    bamfile = bamfn,
    fafile = raer_example("mouse_tiny.fasta"),
    bedfile = raer_example("5k_neuron_sites.bed.gz"),
    cell_barcodes = cbs[1:15],
    verbose = FALSE,
    filterParam = fp
  )
  expect_equal(ncol(se), 15)
  expect_true(all(startsWith(colnames(se), "AAA")))
  expect_true(all(names(assays(se)) == c("nA", "nG")))

  se <- sc_editing(
    bamfile = bamfn,
    fafile = raer_example("mouse_tiny.fasta"),
    bedfile = raer_example("5k_neuron_sites.bed.gz"),
    cell_barcodes = cb_lst,
    verbose = FALSE,
    umi_tag = NULL,
    filterParam = fp
  )
  expect_equal(ncol(se), 5)
  expect_true(all(startsWith(colnames(se), "cluster")))
})


test_that("umi deduplication works", {
  bamfn <- raer_example("5k_neuron_mouse_cbsort.bam")
  idxfn <- build_tag_index(bamfn, overwrite = TRUE)
  cbs <- show_tag_index(bamfn)$tag
  sites <- c("2:598-598", "2:617", "2:645")
  tf <- tempfile(fileext = ".bed")
  on.exit(unlink(tf))
  export(GRanges(sites), tf)


  fp <- FilterParam(library_type = "fr-first-strand",
                    min_base_quality = 1,
                    min_mapq = 0)
  se_umi <- sc_editing(
    bamfile = bamfn,
    fafile = raer_example("mouse_tiny.fasta"),
    bedfile = tf,
    cell_barcodes = cbs[1:25],
    umi_tag = "UB",
    verbose = FALSE,
    filterParam = fp
  )
  se_noumi <- sc_editing(
    bamfile = bamfn,
    fafile = raer_example("mouse_tiny.fasta"),
    bedfile = tf,
    cell_barcodes = cbs[1:25],
    umi_tag = NULL,
    verbose = FALSE,
    filterParam = fp
  )

  av <- as.matrix(assays(se_umi)$nA)
  nz_vals <- av[which(av > 0, arr.ind = TRUE)]
  expect_true(all(nz_vals == 1L))

  av <- as.matrix(assays(se_noumi)$nA)
  nz_vals <- av[which(av > 0, arr.ind = TRUE)]
  expect_true(all(unique(nz_vals) == c(1L, 4L)))

})
