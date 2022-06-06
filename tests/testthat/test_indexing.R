context("indexing")
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)

bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
idx <- indexBed(bedfn)

cbbam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
ubbam_fn <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_ubsort.bam", package = "raer")

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

test_that("tag index works", {
  idx_fn <- paste0(cbbam_fn, ".bri")
  unlink(idx_fn)
  idx_fn <- build_tag_index(cbbam_fn)
  expect_true(file.exists(idx_fn))
  unlink(idx_fn)

  expect_message(build_tag_index(cbbam_fn, n_records_to_check = 0))
  # doesn't exist in first n_records_to_check
  expect_error(build_tag_index(cbbam_fn, tag = "AB"))
  # wrong type
  expect_error(build_tag_index(cbbam_fn, tag = "AS"))
  unlink(idx_fn)
})

test_that("CB tag retrival works", {
    tag_val <- "CB"
    alns <- scanBam(cbbam_fn, param = ScanBamParam(tag = tag_val))
    cbs <- unique(alns[[1]]$tag$CB)

    n_cbs <- sample(seq_along(cbs), 1, replace = FALSE)
    query_cbs <- cbs[seq_len(n_cbs)]
    # idx doesn't exist
    unlink(paste0(cbbam_fn, ".bri"))
    expect_error(get_cell_bam(cbbam_fn, barcodes = cbs))

    idx_fn <- build_tag_index(cbbam_fn, tag = tag_val)
    bam_out <- get_cell_bam(cbbam_fn, barcodes = query_cbs)
    alns <- scanBam(bam_out, param = ScanBamParam(tag = tag_val))

    bam_cbs <- alns[[1]]$tag[[tag_val]]
    expect_equal(length(setdiff(query_cbs, bam_cbs)), 0L,
                 info = paste0("failed with using ", n_cbs, " cbs"))
    expect_true(length(unique(bam_cbs)) == length(query_cbs))

    # output is coordinate sorted
    expect_equal(scanBamHeader(bam_out)[[1]]$text$`@HD`[2] ,
                "SO:coordinate")
    # coordinate sorted bam throws error when indexing
    expect_error(build_tag_index(bam_out))

    unlink(c(bam_out, paste0(bam_out, ".bai")))
    unlink(paste0(cbbam_fn, ".bri"))
})

test_that("other tag (UB) retrival works", {
  tag_val <- "UB"
  alns <- scanBam(ubbam_fn, param = ScanBamParam(tag = tag_val))
  cbs <- unique(alns[[1]]$tag$UB)

  n_cbs <- sample(seq_along(cbs), 1, replace = FALSE)
  query_cbs <- cbs[seq_len(n_cbs)]
  # idx doesn't exist
  unlink(paste0(ubbam_fn, ".bri"))
  expect_error(get_cell_bam(ubbam_fn, barcodes = cbs))

  idx_fn <- build_tag_index(ubbam_fn, tag = tag_val)
  bam_out <- get_cell_bam(ubbam_fn, barcodes = query_cbs)
  alns <- scanBam(bam_out, param = ScanBamParam(tag = tag_val))

  bam_cbs <- alns[[1]]$tag[[tag_val]]
  expect_equal(length(setdiff(query_cbs, bam_cbs)), 0L,
               info = paste0("failed with using ", n_cbs, " cbs"))
  expect_true(length(unique(bam_cbs)) == length(query_cbs))

  unlink(c(bam_out, paste0(bam_out, ".bai")))
  unlink(paste0(ubbam_fn, ".bri"))
})

test_that("invalid input is caught", {
  expect_error(build_tag_index("hello.bam"))
  tmpfn <- tempfile()
  writeLines("hello", tmpfn)
  expect_error(build_tag_index(tmpfn))
  unlink(tmpfn)

  idx_fn <- build_tag_index(cbbam_fn, tag = "CB")
  expect_error(get_cell_bam(cbbam_fn, barcodes = NA))
  expect_error(get_cell_bam(cbbam_fn, barcodes = 1))
  expect_error(get_cell_bam(cbbam_fn, barcodes = character()))

  # barcodes not in bam produce empty bam
  tmp_bam <- get_cell_bam(cbbam_fn, barcodes = "hello")
  alns <- scanBam(tmp_bam)
  expect_true(all(unlist(lapply(alns[[1]], length)) == 0))

  unlink(c(idx_fn, tmp_bam))
})












