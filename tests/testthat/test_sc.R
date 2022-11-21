context("sc_editing")

library(BiocParallel)

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
    filterParam = fp
  )
  expect_equal(ncol(se), 5)
  expect_true(all(startsWith(colnames(se), "cluster")))
})
