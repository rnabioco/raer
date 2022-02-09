context("read_bam")

bam_file <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")

test_that("reading bam works", {
  res <- bam_to_df(bam_file)
  expect_equal(nrow(res), 256)
})


test_that("regional query works", {
  res <- bam_to_df(bam_file, region = "DHFR:11-30")
  expect_equal(nrow(res), 19)
})

test_that("malformed query throws error", {
  expect_error(bam_to_df(bam_file, region = "foo"))
})

test_that("tag extraction works", {
  res <- bam_to_df(bam_file, tags = c("NH:i", "AS:i"))
  cnames <- colnames(res)
  expect_true(all(c("NH", "AS") %in% cnames))
})

test_that("region and tag extraction work", {
  res <- bam_to_df(bam_file, region = "DHFR:11-30",
            tags = c("AS:i"))
  cnames <- colnames(res)
  expect_equal(nrow(res), 19)
  expect_true("AS" %in% cnames))
})

test_that("missing tags are filled with empty strings", {
  res <- bam_to_df(bam_file, tags = "CB:Z")
  expect_true(res$CB[1] == "")
})

test_that("malformed type string throws error", {
  expect_error(bam_to_df(bam_file, tags = "CB:X"))
})

test_that("incorrect types are filled with empty strings", {
  expect_error(all(bam_to_df(bam_file, tags = "NH:Z")$NH == ""))
})
