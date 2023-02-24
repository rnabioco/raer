library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)
bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.md.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")

sites <- import(bedfn)
tmp_dir <- tempdir("plps/tests")
on.exit(unlink(tmp_dir, recursive = TRUE))

out_fns <- pileup_sites(bamfn, fafn, sites,
                        return_data = FALSE,
                        outfile_prefix = tmp_dir)
out_short_fns <- pileup_sites(bamfn, fafn, region = "SSR3:203-205",
                              return_data = FALSE,
                              outfile_prefix = paste0(tmp_dir, "_short"))
plp <- read_pileup(out_fns[2])
plp2 <- plp
plp_short <- read_pileup(out_short_fns[2])

test_that("creating a RangedSummarizedExperiment works", {
  # One sample
  se_obj <- merge_pileups(plp)
  expect_equal(length(assays(se_obj)), 7)
  se_obj <- merge_pileups(plp, sample_names = "sample_1")
  expect_equal(length(assays(se_obj)), 7)

  # Two samples
  expect_error(merge_pileups(plp, sample_names = c("sample_1", "sample_2")))
  se_obj <- merge_pileups(list(plp, plp2))
  expect_equal(length(assays(se_obj)), 7)
  expect_equal(ncol(assays(se_obj)[[1]]), 2)
  expect_error(merge_pileups(list(plp, plp2), c("sample_1")))

  # Three samples, different lengths
  se_obj <- merge_pileups(list(plp, plp_short, plp2))
  expect_equal(length(assays(se_obj)), 7)
  expect_equal(ncol(assays(se_obj)[[1]]), 3)
  expect_equal(nrow(assays(se_obj)[[1]]), 182)
})


