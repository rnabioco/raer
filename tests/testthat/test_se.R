context("create_se")

library(GenomicRanges)
library(SummarizedExperiment)

bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
plp <- get_pileup(bamfn, fafn, bedfn)
plp2 <- get_pileup(bamfn, fafn, bedfn)
plp3 <- get_pileup(bamfn, fafn, bedfn)
plp_short <- get_pileup(bamfn, fafn, region = "SSR3:203-205")

test_that("creating a RangedSummarizedExperiment works", {
  # One sample
  se_obj <- create_se(plp)
  expect_equal(length(assays(se_obj)), 7)
  se_obj_2 <- create_se(plp, sample_names = "sample_1")
  expect_equal(length(assays(se_obj)), 7)

  # Two samples
  expect_error(create_se(plp, sample_names =c("sample_1", "sample_2")))
  se_obj_3 <- create_se(list(plp, plp2))
  expect_equal(length(assays(se_obj_3)), 7)
  expect_equal(ncol(assays(se_obj_3)[[1]]), 2)
  expect_error(create_se(list(plp, plp2), c("sample_1")))
  se_obj_4 <- create_se(list(plp, plp2), sample_names = c("sample_1", "sample_2"))
  expect_equal(length(assays(se_obj_4)), 7)
  expect_equal(ncol(assays(se_obj_4)[[1]]), 2)

  # Three samples
  se_obj5 <- create_se(list(plp, plp2, plp3),
                       sample_names = c("sample_1", "sample_2", "sample_3"))
  expect_equal(length(assays(se_obj5)), 7)
  expect_equal(ncol(assays(se_obj5)[[1]]), 3)

  # Three samples, different lengths
  se_obj_6 <- create_se(list(plp, plp_short, plp2))
  expect_equal(length(assays(se_obj5)), 7)
  expect_equal(ncol(assays(se_obj5)[[1]]), 3)
  expect_equal(nrow(assays(se_obj5)[[1]]), 182)
})


test_that("creating a sparse RangedSummarizedExperiment works", {
  # non-numeric cols throw warning
  expect_warning(sres <- create_se(plp, sparse = TRUE))
  expect_true(is(assay(sres, "nA"), "sparseMatrix"))
  expect_false(is(assay(sres, "Var"), "sparseMatrix"))
  dres <- create_se(plp, sparse = FALSE)

  assays(sres) <- lapply(assays(sres), function(x){
    if(is(x, "sparseMatrix")){
      x <- as.matrix(x)
      mode(x) <- 'integer'
    }
    x})

  expect_true(identical(assays(dres), assays(sres)))

  plps <- list(plp, plp_short, plp2)
  # sparse coerces all missing values to 0
  sres <- create_se(plps, assay_cols = c("nA", "nG"), sparse = TRUE)
  dres <- create_se(plps, assay_cols = c("nA", "nG"),
                    sparse = FALSE,
                    fill_na = 0L)

  expect_true(all(names(assays(sres)) == c("nA", "nG")))

  assays(sres) <- lapply(assays(sres), function(x){
    if(is(x, "sparseMatrix")){
      x <- as.matrix(x)
      mode(x) <- 'integer'
    }
    x})
  expect_true(identical(assays(dres), assays(sres)))
})
