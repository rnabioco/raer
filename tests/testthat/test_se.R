context("create_se")

library(GenomicRanges)
library(SummarizedExperiment)

bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
res <- get_pileup(bamfn, fafn, bedfn)
res2 <- get_pileup(bamfn, fafn, bedfn)
res3 <- get_pileup(bamfn, fafn, bedfn)
res_short <- get_pileup(bamfn, fafn, region = "SSR3:203-205")

test_that("creating a RangedSummarizedExperiment works", {
  # One sample
  se_obj <- create_se(res)
  expect_equal(length(assays(se_obj)), 7)
  se_obj_2 <- create_se(res, sample_names = "sample_1")
  expect_equal(length(assays(se_obj)), 7)

  # Two samples
  expect_error(create_se(res, sample_names =c("sample_1", "sample_2")))
  se_obj_3 <- create_se(list(res, res2))
  expect_equal(length(assays(se_obj_3)), 7)
  expect_equal(ncol(assays(se_obj_3)[[1]]), 2)
  expect_error(create_se(list(res, res2), c("sample_1")))
  se_obj_4 <- create_se(list(res, res2), sample_names = c("sample_1", "sample_2"))
  expect_equal(length(assays(se_obj_4)), 7)
  expect_equal(ncol(assays(se_obj_4)[[1]]), 2)

  # Three samples
  se_obj5 <- create_se(list(res, res2, res3),
                       sample_names = c("sample_1", "sample_2", "sample_3"))
  expect_equal(length(assays(se_obj5)), 7)
  expect_equal(ncol(assays(se_obj5)[[1]]), 3)

  # Three samples, different lengths
  se_obj_6 <- create_se(list(res, res_short, res2))
  expect_equal(length(assays(se_obj5)), 7)
  expect_equal(ncol(assays(se_obj5)[[1]]), 3)
  expect_equal(nrow(assays(se_obj5)[[1]]), 182)
})

