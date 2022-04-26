context("indexing")
library(GenomicRanges)

bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
idx <- indexBed(bedfn)

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

