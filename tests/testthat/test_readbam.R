context("read_bam")

bam_file <- system.file("extdata", "small_sorted.bam", package = "kentr")

test_that("reading bam works", {
  res <- bam_to_df(bam_file)
  expect_equal(nrow(res), 9976)
})


test_that("regional query works", {
  res <- bam_to_df(bam_file, region = "chr1:3e6-3.1e6")
  expect_equal(nrow(res), 16)
})

test_that("malformed query throws error", {
  expect_error(bam_to_df(bam_file, region = "foo"))
})

test_that("tag extraction works", {
  res <- bam_to_df(bam_file, tags = c("XS:A", "AS:i"))
  cnames <- colnames(res)
  expect_true(all(c("XS", "AS") %in% cnames))
})

test_that("region and tag extraction work", {
  res <- bam_to_df(bam_file, region = "chr1:3e6-3.1e6",
            tags = c("XS:A", "AS:i"))
  cnames <- colnames(res)
  expect_equal(nrow(res), 16)
  expect_true(all(c("XS", "AS") %in% cnames))
})

test_that("missing tags are filled with empty strings", {
  res <- bam_to_df(bam_file, tags = "CB:Z")
  expect_true(res$CB[1] == "")
})

test_that("malformed type string throws error", {
  expect_error(bam_to_df(bam_file, tags = "CB:X"))
})
