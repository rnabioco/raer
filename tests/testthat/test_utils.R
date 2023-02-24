
test_that("get_region works", {
  x <- get_region("chr1:1-20")
  expect_true(all(names(x) == c("chrom", "start", "end")))
  expect_true(x$chrom == "chr1" && x$start == 0 && x$end == 20)

  x <- get_region("chr1")
  expect_true(x$start == 0 && x$end == 2^31 - 1)
  expect_error(get_region("214124:!124:1:!24:chr1"))
})


