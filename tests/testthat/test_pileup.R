context("get_pileup")
library(GenomicRanges)

bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.bam", package = "raer")
# bamfn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")


fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")

nts <- c("A", "T", "G", "C")
nt_clmns <- paste0("n", nts)

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
  expect_equal(length(res$Ref), 529)
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

test_that("pileup chrom start end args", {
  res <- get_pileup(
    bamfn, fafn, chrom = "DHFR"
  )

  expect_identical(levels(res@seqnames), "DHFR")
})

test_that("pileup depth lims", {
  unflt <- get_pileup(bamfn, fafn)
  unflt <- as.data.frame(unflt)
  unflt$seqnames <- as.character(unflt$seqnames)

  rsums <- unflt[nt_clmns]
  rsums <- rowSums(rsums) >= 30

  expected_df <- unflt[rsums, ]
  rownames(expected_df) <- NULL

  res <- get_pileup(bamfn, fafn, min_reads = 30)
  res <- as.data.frame(res, row.names = NULL)
  res$seqnames <- as.character(res$seqnames)

  expect_equal(expected_df, res)

  # How does max_depth filter?
  # res <- get_pileup(bamfn, fafn, max_depth = 10)
  # res <- as.data.frame(res)
  # rsums <- res[nt_clmns]
  # rsums <- rowSums(rsums)
})





check_nRef_calc <- function(input, nts_in = nts) {
  res <- as.data.frame(input)

  for (nt in nts_in) {
    clmn <- paste0("n", nt)
    dat  <- res[res$Ref == nt, ]

    expect_identical(dat[, clmn], dat$nRef)

    other_clmns <- nts_in[nts_in != nt]
    other_clmns <- paste0("n", other_clmns)
    var_sums    <- rowSums(dat[, other_clmns])

    expect_identical(as.integer(var_sums), as.integer(dat$nVar))
  }
}

test_that("pileup check nRef and nVar", {

  # strds <- c("unstranded", "fr-first-strand", "fr-second-strand")
  strds <- c("fr-first-strand", "fr-second-strand")

  for (strd in strds) {
    res <- get_pileup(bamfn, fafn, library_type = strd)

    check_nRef_calc(res)
  }
})





# TEST CALCULATIONS
res1 <- get_pileup(
  bamfn, fafn, min_reads = 0,
  min_base_qual = 0,
  library_type = "fr-first-strand"
) %>%
  as.data.frame()

res2 <- get_pileup(
  bamfn, fafn, min_reads = 0,
  min_base_qual = 0,
  library_type = "fr-second-strand"
) %>%
  as.data.frame()

# nt <- c("A", "T")
#
# clmn  <- paste0("n", nt)
# dat   <- res[res$Ref %in% nt, ]
# rsums <- rowSums(dat[, clmn])
#
# expect_equal(unname(rsums), dat$nRef)
#
# other_clmns <- nts_in[nts_in != nt]
# other_clmns <- paste0("n", other_clmns)
# var_sums    <- rowSums(dat[, other_clmns])
#
# expect_identical(as.integer(var_sums), as.integer(dat$nVar))





