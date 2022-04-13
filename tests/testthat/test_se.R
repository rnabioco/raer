context("create_se")

library(GenomicRanges)
library(SummarizedExperiment)

bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
res <- get_pileup(bamfn, fafn, bedfn)
res2 <- get_pileup(bamfn, fafn, bedfn)
res3 <- get_pileup(bamfn, fafn, bedfn)
res_short <- get_pileup(bamfn, fafn, region = "SSR3:203-205")

test_that("creating a RangedSummarizedExperiment works", {
  # One sample
  se_obj <- create_se(res)
  expect_equal(length(assays(se_obj)), 9)
  se_obj_2 <- create_se(res, "sample_1")
  expect_equal(length(assays(se_obj)), 9)

  # Two samples
  expect_error(create_se(res, c("sample_1", "sample_2")))
  se_obj_3 <- create_se(list(res, res2))
  expect_equal(length(assays(se_obj_3)), 9)
  expect_equal(ncol(assays(se_obj_3)$Ref), 2)
  expect_error(create_se(list(res, res2), c("sample_1")))
  se_obj_4 <- create_se(list(res, res2), c("sample_1", "sample_2"))
  expect_equal(length(assays(se_obj_4)), 9)
  expect_equal(ncol(assays(se_obj_4)$Ref), 2)

  # Three samples
  se_obj5 <- create_se(list(res, res2, res3),
                       c("sample_1", "sample_2", "sample_3"))
  expect_equal(length(assays(se_obj5)), 9)
  expect_equal(ncol(assays(se_obj5)$Ref), 3)

  # Three samples, different lengths
  se_obj_6 <- create_se(list(res, res_short, res2))
  expect_equal(length(assays(se_obj5)), 9)
  expect_equal(ncol(assays(se_obj5)$Ref), 3)
  expect_equal(nrow(assays(se_obj5)$Ref), 182)
})

