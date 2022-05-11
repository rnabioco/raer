context("indexing")
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)

bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
idx <- indexBed(bedfn)

cbbam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")

test_that("bedindex works", {
  expect_true(is(idx, "BedFile"))
  expect_false(is_null_extptr(idx$.extptr))
})

test_that("bedindex can be used with pileup", {
  res <- get_pileup(bamfn, fafn, bedidx = idx)
  res2 <- get_pileup(bamfn, fafn, bedfile = bedfn)
  expect_true(identical(res, res2))
})

test_that("close cleans up index", {
  idx <- close(idx)
  expect_false(idx$open)
  expect_error(get_pileup(bamfn, fafn,  bedidx = idx))
})

# test_that("tag index works", {
#   idx_fn <- paste0(cbbam_fn, ".bri")
#   unlink(idx_fn)
#   idx_fn <- build_tag_index(cbbam_fn)
#   idx_fn <- build_tag_index(cbbam_fn)
#   expect_true(file.exists(idx_fn))
#   unlink(idx_fn)
#
# })
#
# test_that("CB retrival works", {
#
#     alns <- readGAlignments(cbbam_fn, param = ScanBamParam(tag = "CB"))
#     cbs <- unique(mcols(alns)$CB)
#
#     n_cbs <- sample(seq_along(cbs), 1, replace = FALSE)
#     query_cbs <- cbs[seq_len(n_cbs)]
#     # idx doesn't exist
#     unlink(paste0(cbbam_fn, ".bri"))
#     expect_error(get_cell_bam(cbbam_fn, barcodes = cbs))
#
#     idx_fn <- build_tag_index(cbbam_fn)
#     bam_out <- get_cell_bam(cbbam_fn, barcodes = query_cbs)
#     bam_out <- get_cell_bam(cbbam_fn, barcodes = query_cbs)
#     alns <- readGAlignments(bam_out, param = ScanBamParam(tag = "CB"))
#
#     bam_cbs <- mcols(alns)$CB
#     expect_equal(length(setdiff(query_cbs, bam_cbs)), 0L,
#                  info = paste0("failed with using ", n_cbs, " cbs"))
#     expect_true(length(unique(bam_cbs)) == length(query_cbs))
#
#     # output is coordinate sorted
#     expect_equal(scanBamHeader(bam_out)[[1]]$text$`@HD`[2] ,
#                 "SO:coordinate")
#
#     unlink(c(bam_out, paste0(bam_out, ".bai")))
# })

