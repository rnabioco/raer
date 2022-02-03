context("get_pileup")
library(GenomicRanges)

bamfn <- system.file("extdata", "SRR1258218_chr21.bam", package = "ullr")
fafn <- system.file("extdata", "chr21.fa", package = "ullr")
bedfn <- system.file("extdata", "chr21_regions.bed", package = "ullr")
res <- get_pileup(bamfn, fafn, bedfn)

test_that("pileup works", {
  expect_equal(length(res$Ref), 5)
})

test_that("pileup regional query works", {

  res <- get_pileup(bamfn, fafn, region = "chr21:46664718-46664718")
  expect_equal(length(res$Ref), 1)
  expect_equal(start(res), 46664718)
  expect_equal(end(res), 46664719)

  # chr1 does not exist
  res <- get_pileup(bamfn, fafn, region = "chr1")
  expect_null(res)

  res <- get_pileup(bamfn, fafn, bedfile = NULL, chrom = "chr21")
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
  res_first <- get_pileup(bamfn, fafn, bedfn, library_type = "fr-first-strand")
  res_second <- get_pileup(bamfn, fafn, bedfn, library_type = "fr-second-strand")
  res_un <- get_pileup(bamfn, fafn, bedfn, library_type = "unstranded")
  expect_true(all(strand(res_first) == "+"))
  expect_true(all(strand(res_second) == "-"))
  expect_true(all(strand(res_un) == "-"))
  expect_error(get_pileup(bamfn, fafn, bedfn, library_type = "unknown-string"))
})