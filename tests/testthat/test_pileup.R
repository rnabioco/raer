context("get_pileup")
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)

bamfn <- system.file("extdata", "SRR5564269_Aligned.sortedByCoord.out.bam", package = "raer")
bam2fn <- system.file("extdata", "SRR5564277_Aligned.sortedByCoord.out.bam", package = "raer")
fafn <- system.file("extdata", "human.fasta", package = "raer")
bedfn <- system.file("extdata", "regions.bed", package = "raer")

nts <- c("A", "T", "G", "C")
nt_clmns <- paste0("n", nts)

res <- get_pileup(bamfn, fafn, bedfn)

test_that("pileup works", {
  expect_equal(length(res$Ref), 182)
  expect_equal(ncol(as.data.frame(res)), 14)
})

test_that("filtering for variants in pileup works", {
  vars <- res[res$Var != "-"]
  seqlevels(vars) <- seqlevelsInUse(vars)
  res_all_vars <- get_pileup(bamfn, fafn, bedfn, only_keep_variants = TRUE)
  expect_equal(vars, res_all_vars)
})

test_that("2-bam pileup works", {
  count_cols <- c("nRef", "nVar", "nA", "nT", "nC", "nG", "nN")
  res <- get_pileup(c(bamfn, bamfn), fafn, bedfn,
                    library_type = c("fr-first-strand",
                                     "fr-first-strand"))

  expect_equal(length(res$Ref), 182)
  expect_equal(ncol(as.data.frame(res)), 21)
  b1_vals <- mcols(res)[paste0(count_cols, "_1")]
  b2_vals <- mcols(res)[paste0(count_cols, "_2")]
  colnames(b2_vals) <- colnames(b1_vals)
  expect_true(identical(b1_vals, b2_vals))

  # should default to fr-first-strand and unstranded
  res <- get_pileup(c(bamfn, bamfn), fafn, bedfn)

  expect_equal(length(res$Ref), 182)
  expect_equal(ncol(as.data.frame(res)), 21)
  b1_vals <- mcols(res)[paste0(count_cols, "_1")]
  b2_vals <- mcols(res)[paste0(count_cols, "_2")]
  colnames(b2_vals) <- colnames(b1_vals)
  expect_false(identical(b1_vals, b2_vals))

  res <- get_pileup(c(bamfn, bam2fn), fafn, bedfn)
  b1_vals <- mcols(res)[paste0(count_cols, "_1")]
  b2_vals <- mcols(res)[paste0(count_cols, "_2")]
  colnames(b2_vals) <- colnames(b1_vals)
  expect_false(identical(b1_vals, b2_vals))
})


test_that("pileup regional query works", {

  res <- get_pileup(bamfn, fafn, region = "SSR3:203-205")
  expect_equal(length(res$Ref), 3)
  expect_equal(start(res), c(203, 204, 205))
  expect_equal(end(res), c(203, 204, 205))

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
  res <- get_pileup(bamfn, fafn, library_type = "fr-first-strand")
  expect_true(table(strand(res))["+"] == 619)
  expect_true(table(strand(res))["-"] == 1047)
  expect_true(all(strand(res[seqnames(res) == "DHFR"]) == "-"))
  expect_true(all(strand(res[seqnames(res) == "SPCS3"]) == "+"))

  res <- get_pileup(bamfn, fafn, library_type = "fr-second-strand")
  expect_true(table(strand(res))["+"] == 1047)
  expect_true(table(strand(res))["-"] == 619)
  expect_true(all(strand(res[seqnames(res) == "DHFR"]) == "+"))
  expect_true(all(strand(res[seqnames(res) == "SPCS3"]) == "-"))

  res <- get_pileup(bamfn, fafn, library_type = "unstranded")
  expect_true(all(strand(res) == "+"))

  expect_error(get_pileup(bamfn, fafn, bedfn, library_type = "unknown-string"))
})

test_that("pileup chrom start end args", {
  res <- get_pileup(
    bamfn, fafn, chrom = "DHFR"
  )

  expect_identical(seqlevels(res), "DHFR")
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

bout <- tempfile(fileext = ".bam")

test_that("pileup check nh tag filter", {
  bout <- filterBam(bamfn,
                    param = ScanBamParam(tagFilter = list("NH"= 1)),
                    destination = bout,
                    indexDestination = TRUE)
  a <- get_pileup(bamfn, fafn, n_align = 1, n_align_tag = "NH")
  b <- get_pileup(bout, fafn)
  expect_true(identical(a, b))
})

test_that("pileup check mapq filter", {
  bout <- filterBam(bamfn,
                    param = ScanBamParam(mapqFilter = 255),
                    destination = bout,
                    indexDestination = TRUE)
  a <- get_pileup(bamfn, fafn, min_mapq = 255)
  b <- get_pileup(bout, fafn)
  expect_true(identical(a, b))
})

test_that("pileup check mapq filter", {
  bout <- filterBam(bamfn,
                    param = ScanBamParam(mapqFilter = 255),
                    destination = bout,
                    indexDestination = TRUE)
  a <- get_pileup(bamfn, fafn, min_mapq = 255)
  b <- get_pileup(bout, fafn)
  expect_true(identical(a, b))
})


test_that("pileup check flag filtering", {
  bout <- filterBam(bamfn,
                    param = ScanBamParam(flag = scanBamFlag(isMinusStrand = F)),
                    destination = bout,
                    indexDestination = TRUE)
  a <- get_pileup(bamfn, fafn, bedfn, bam_flags = scanBamFlag(isMinusStrand = F))
  b <- get_pileup(bout, fafn, bedfn)
  expect_true(identical(a, b))

  bout <- filterBam(bamfn,
                    param = ScanBamParam(flag = scanBamFlag(isMinusStrand = T)),
                    destination = bout,
                    indexDestination = TRUE)
  a <- get_pileup(bamfn, fafn, bedfn, bam_flags = scanBamFlag(isMinusStrand = T))
  b <- get_pileup(bout, fafn, bedfn)
  expect_true(identical(a, b))
})

fa <- scanFa(fafn)
hp_matches <- vmatchPattern(strrep("A", 6), fa) %>% as(., "GRanges")

test_that("pileup check homopolymer filter", {
  a <- get_pileup(bamfn, fafn, event_filters = c(0, 0, 0, 6))
  b <- get_pileup(bamfn, fafn, event_filters = c(0, 0, 0, 0))
  expect_false(identical(a, b))

  expect_equal(length(queryHits(findOverlaps(hp_matches, a))), 0)
  expect_equal(length(queryHits(findOverlaps(hp_matches, b))), 120)
})

test_that("filtering for splicing events works",{
  # get sites near splicing events
  splices <- GRanges(coverage(junctions(GenomicAlignments::readGAlignmentPairs(bamfn))))
  splices <- splices[splices$score > 0]
  plp <- get_pileup(bamfn, fafn,
                  event_filters = c(0, 0, 0, 0))
  plp <- plp[plp$Var != "-"]

  # 3 sites, 1 is from all non-spliced reds, 2 sites are all from spliced reads,
  # (SPCS3  227, 347, 348)
  sites_near_splices <- plp[queryHits(findOverlaps(plp, splices, maxgap = 4))]

  plp_b <- get_pileup(bamfn, fafn, event_filters = c(0, 5, 0, 0))
  plp_b <- plp_b[plp_b$Var != "-"]
  bsites_near_splices <- plp_b[queryHits(findOverlaps(plp_b,
                                                      splices,
                                                      maxgap = 4))]
  expect_equal(length(sites_near_splices$Ref), 3)
  expect_equal(length(bsites_near_splices$Ref), 1)
})


test_that("filtering for indel events works",{
  # get sites near splicing events
  reads <- GenomicAlignments::readGAlignments(bamfn, use.names = T)
  cig_ops <- cigarRangesAlongReferenceSpace(cigar(reads), ops = "D", pos = start(reads))
  refs <- seqnames(reads)[elementNROWS(cig_ops) > 0]
  cig_ops <- cig_ops[elementNROWS(cig_ops) > 0]
  names(cig_ops) <- refs
  cig_pos <- unique(GRanges(cig_ops))

  plp <- get_pileup(bamfn, fafn,
                    event_filters = c(0, 0, 0, 0),
                    min_base_qual = 1)
  plp <- plp[plp$Var != "-"]

  # SSR3 387 variant has 1 read with variant with indel 4 bp away
  # SRR3 388 variant has 1 read with variant 3bp away from indel
  sites_near_indels <- plp[queryHits(findOverlaps(plp, cig_pos))]

  plp_b <- get_pileup(bamfn, fafn, event_filters = c(0, 0, 5, 0),
                      min_base_qual = 1)
  plp_b <- plp_b[plp_b$Var != "-"]
  bsites <- plp_b[queryHits(findOverlaps(plp_b,
                                         cig_pos))]

  expect_equal(length(sites_near_indels$Ref), 3)
  # lose SSR3 387 site after filtering
  expect_equal(length(bsites$Ref), 2)
  # lose 1 read from SSR3 388 site after filtering
  expect_equal(sites_near_indels$nVar[2] - bsites$nVar[1], 1)
})

bout_sorted <- tempfile()
bed <- tempfile(fileext = ".bed")
test_that("writing reads with mismatches works", {
  plp <- get_pileup(bamfn, fafn, outbam = bout)
  plp <- plp[plp$Var != "-"]
  # note that this removes single end reads, need to consider how to handle these
  bam <- GenomicAlignments::readGAlignmentPairs(bout, strandMode = 2)
  # check that all reads overlap sites with variants
  overlapping_alignments = bam[unique(queryHits(findOverlaps(bam, plp)))]
  expect_identical(bam, overlapping_alignments)

  # check that rerunning with bam file of variant overlapping reads
  # reproduces # and type of variants
  srted_bam_name <- Rsamtools::sortBam(bout, bout_sorted)
  export.bed(plp, bed)
  replp <- get_pileup(srted_bam_name, fafn, bedfile = bed)

  expect_equal(ranges(plp), ranges(replp))
  expect_equal(plp$nVar, replp$nVar)
  expect_equal(plp$Var, replp$Var)
})


unlink(c(bout, bout_sorted, bed))
