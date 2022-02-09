context("get_pileup")
library(GenomicRanges)

bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")
res <- get_pileup(bamfn, fafn, bedfn)

test_that("pileup works", {
  expect_equal(length(res$Ref), 182)
})

test_that("pileup regional query works", {

  res <- get_pileup(bamfn, fafn, region = "SSR3:203-205")
  expect_equal(length(res$Ref), 3)
  expect_equal(start(res), c(203, 204, 205))
  expect_equal(end(res), c(204, 205, 206))

  # chr1 does not exist
  expect_error(get_pileup(bamfn, fafn, region = "chr1"))

  res <- get_pileup(bamfn, fafn, bedfile = NULL, chrom = "SSR3")
  expect_true(length(res$Ref), 529)
})

test_that("incorrect regional query is caught", {
  expect_error(get_pileup(bamfn, fafn, region = "chrHello"))
})

test_that("missing files are caught", {
  expect_error(get_pileup(bamfile = "hello.bam", fafn, bedfn))
  expect_error(get_pileup(bamfn, fafile = "hello.fasta", bedfn))
  expect_error(get_pileup(bamfn, fafn, bedfile = "hello.bed"))
})

test_that("library types are respected", {
  res_first <- get_pileup(bamfn, fafn, library_type = "fr-first-strand")
  res_second <- get_pileup(bamfn, fafn, library_type = "fr-second-strand")
  res_un <- get_pileup(bamfn, fafn, library_type = "unstranded")
  # ... add tests for strands here.
  expect_error(get_pileup(bamfn, fafn, bedfn, library_type = "unknown-string"))
})