library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(BiocParallel)
library(stringr)
library(rtracklayer)

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
bedfn <- raer_example("regions.bed")
sites <- import(bedfn)

res <- pileup_sites(bamfn, fafn, sites)

test_that("pileup works", {
  expect_equal(length(rowData(res)$REF), 182)
  expect_equal(length(assays(res)), 7)
})

test_that("filtering for variants in pileup works", {
  vars <- res[assay(res, "ALT")[, 1] != "-"]
  res_all_vars <- pileup_sites(bamfn, fafn, sites,
    param = FilterParam(only_keep_variants = TRUE)
  )
  expect_equal(nrow(res_all_vars), 2)
  expect_equal(vars, res_all_vars)

  res_2_vars <- pileup_sites(bamfn, fafn, sites,
                             param = FilterParam(min_variant_reads = 2))
  expect_equal(nrow(res_2_vars), 1)
  expect_equal(vars[assay(vars, "nAlt")[, 1] >= 2, ],
               res_2_vars)

})

test_that("n-bam pileup works", {
  res <- pileup_sites(c(bamfn, bamfn), fafn, sites,
    param = FilterParam(library_type = c(
      "fr-first-strand",
      "fr-first-strand"
    ))
  )
  expect_true(identical(res[, 1], res[, 2]))
  expect_equal(dim(res), c(182, 2))
  expect_equal(length(assays(res)), 7)

  res <- pileup_sites(c(bamfn, bam2fn), fafn, sites)
  expect_equal(dim(res), c(182, 2))
  expect_false(identical(res[, 1], res[, 2]))

  res <- pileup_sites(c(bamfn, bamfn), fafn, sites,
    param = FilterParam(library_type = c(
      "fr-first-strand",
      "genomic-unstranded"
    ))
  )
  expect_equal(dim(res), c(214, 2))
  no_nas <- all(vapply(assays(res), function(x) all(!is.na(x)), logical(1)))
  expect_true(no_nas)

  res <- pileup_sites(c(bamfn, bam2fn), fafn, sites,
    param = FilterParam(
      library_type = "fr-first-strand",
      only_keep_variants = TRUE
    ),
  )
  # same sites reported in both files
  expect_true(all(start(res[, 1]) == start(res[, 2])))
  # sites are variant in at least 1 file
  expect_true(all(rowSums(assay(res, "ALT") != "-") > 0))

  res <- pileup_sites(rep(bamfn, 4), fafn, sites,
    param = FilterParam(library_type = "fr-first-strand")
  )
  expect_true(ncol(res) == 4)

  # all sites in first bam are variant
  res <- pileup_sites(c(bamfn, bam2fn), fafn, sites,
    param = FilterParam(
      library_type = "fr-first-strand",
      only_keep_variants = c(TRUE, FALSE)
    )
  )
  expect_true(all(assay(res[,1], "ALT") != "-"))

  # all sites in second bam are variant
  res <- pileup_sites(c(bamfn, bam2fn), fafn, sites,
    param = FilterParam(
      library_type = "fr-first-strand",
      only_keep_variants = c(FALSE, TRUE)
    )
  )
  expect_true(all(assay(res[,2], "ALT") != "-"))
})

test_that("pileup regional query works", {
  res <- pileup_sites(bamfn, fafn, region = "SSR3:203-205")
  expect_equal(nrow(res), 3)
  expect_equal(start(res), c(203, 204, 205))
  expect_equal(end(res), c(203, 204, 205))

  # chr1 does not exist
  expect_error(suppressWarnings(pileup_sites(bamfn, fafn, region = "chr1")))

  res <- pileup_sites(bamfn, fafn, chrom = "SSR3")
  expect_equal(nrow(res), 529)
})

test_that("incorrect regional query is caught", {
  # will produce warning that chrHello is not in bamfile and an error
  expect_error(suppressWarnings(pileup_sites(bamfn, fafn, region = "chrHello")))
})

test_that("missing files are caught", {
  expect_error(pileup_sites(bamfile = "hello.bam", fafn))
  expect_error(pileup_sites(bamfn, fafile = "hello.fasta"))
})

test_that("library types are respected", {
  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(library_type = "fr-first-strand")
  )
  expect_true(table(strand(res))["+"] == 619)
  expect_true(table(strand(res))["-"] == 1047)
  expect_true(all(strand(res[seqnames(res) == "DHFR"]) == "-"))
  expect_true(all(strand(res[seqnames(res) == "SPCS3"]) == "+"))

  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(library_type = "fr-second-strand")
  )
  expect_true(table(strand(res))["+"] == 1047)
  expect_true(table(strand(res))["-"] == 619)
  expect_true(all(strand(res[seqnames(res) == "DHFR"]) == "+"))
  expect_true(all(strand(res[seqnames(res) == "SPCS3"]) == "-"))

  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(library_type = "genomic-unstranded")
  )
  expect_true(all(strand(res) == "+"))

  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(library_type = "unstranded")
  )
  expect_true(all(strand(res) %in% c("+", "-")))

  expect_error(pileup_sites(bamfn, fafn, sites,
    param = FilterParam(library_type = "unknown-string")
  ))
})

test_that("pileup chrom start end args", {
  res <- pileup_sites(
    bamfn, fafn,
    chrom = "DHFR"
  )
  expect_true(all(seqnames(res) == "DHFR"))
})

nts <- c("A", "T", "G", "C")
nt_clmns <- paste0("n", nts)

test_that("pileup depth lims", {
  unflt <- pileup_sites(bamfn, fafn)
  base_cnts <- do.call(cbind, assays(unflt)[nt_clmns])
  flt_cnts <- base_cnts[rowSums(base_cnts) >= 30, ]

  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(min_depth = 30)
  )
  base_cnts <- do.call(cbind, assays(res)[nt_clmns])
  expect_equal(base_cnts, base_cnts)
})

check_nRef_calc <- function(input, nts_in = nts) {
  res <- as.data.frame(do.call(cbind, assays(input)))
  colnames(res) <- names(assays(input))
  for (nt in nts_in) {
    clmn <- paste0("n", nt)
    dat <- res[res$REF == nt, ]

    expect_identical(dat[, clmn], dat$nRef)

    other_clmns <- nts_in[nts_in != nt]
    other_clmns <- paste0("n", other_clmns)
    var_sums <- rowSums(dat[, other_clmns])

    expect_identical(as.integer(var_sums), as.integer(dat$nAlt))
  }
}

test_that("pileup check nRef and nAlt", {
  strds <- c("fr-first-strand", "fr-second-strand",
             "unstranded", "genomic-unstranded")

  for (strd in strds) {
    res <- pileup_sites(bamfn, fafn,
                      param = FilterParam(library_type = strd))
    check_nRef_calc(res)
  }
})

bout <- tempfile(fileext = ".bam")

test_that("pileup check mapq filter", {
  bout <- filterBam(bamfn,
    param = ScanBamParam(mapqFilter = 255),
    destination = bout,
    indexDestination = TRUE
  )
  # enforce same colname in rse returned
  names(bout) <- basename(bamfn)
  a <- pileup_sites(bamfn, fafn, param = FilterParam(min_mapq = 255))
  b <- pileup_sites(bout, fafn)
  expect_true(identical(a, b))
})

test_that("pileup check flag filtering", {
  bout <- filterBam(bamfn,
    param = ScanBamParam(flag = scanBamFlag(isMinusStrand = F)),
    destination = bout,
    indexDestination = TRUE
  )
  names(bout) <- basename(bamfn)
  fp <- FilterParam(bam_flags = scanBamFlag(isMinusStrand = F))
  a <- pileup_sites(bamfn, fafn, sites, param = fp)
  b <- pileup_sites(bout, fafn, sites)
  expect_true(identical(a, b))

  bout <- filterBam(bamfn,
    param = ScanBamParam(flag = scanBamFlag(isMinusStrand = T)),
    destination = bout,
    indexDestination = TRUE
  )
  names(bout) <- basename(bamfn)
  fp <- FilterParam(bam_flags = scanBamFlag(isMinusStrand = T))
  a <- pileup_sites(bamfn, fafn, sites, param = fp)
  b <- pileup_sites(bout, fafn, sites)
  expect_true(identical(a, b))
})



test_that("pileup read trimming filter works", {
  vec <- function(rse, assay, col = 1) unname(assays(rse)[[assay]][, col])

  a <- pileup_sites(bamfn, fafn, region = "SSR3:440-450")
  b <- pileup_sites(bamfn, fafn, param = FilterParam(trim_5p = 6), region = "SSR3:440-450")
  d <- pileup_sites(bamfn, fafn, param = FilterParam(trim_5p = 6, trim_3p = 6), region = "SSR3:440-450")
  e <- pileup_sites(bamfn, fafn, param = FilterParam(trim_3p = 6), region = "SSR3:440-450")
  ex1 <- c(0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L)
  ex2 <- c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)

  expect_equal(vec(a, "nRef") -  vec(b, "nRef"), ex1)
  expect_equal(vec(b, "nRef") - vec(d, "nRef"), ex2)
  expect_equal(vec(a, "nRef") - vec(e, "nRef"), ex2)

  a <- pileup_sites(bamfn, fafn, region = "SSR3:128-130")
  b <- pileup_sites(bamfn, fafn, param = FilterParam(trim_5p = 6), region = "SSR3:128-130")
  d <- pileup_sites(bamfn, fafn, param = FilterParam(trim_5p = 6, trim_3p = 6), region = "SSR3:128-130")
  e <- pileup_sites(bamfn, fafn, param = FilterParam(trim_3p = 6), region = "SSR3:128-130")
  ex1 <- c(3L, 4L, 2L)
  ex2 <- c(3L, 2L, 0L)
  expect_equal(vec(a, "nRef") - vec(b, "nRef"), ex1)
  expect_equal(vec(b, "nRef") - vec(d, "nRef"), ex2)
  expect_equal(vec(a, "nRef") - vec(e, "nRef"), ex2)
})

test_that("pileup fractional read trimming filter works", {
  vec <- function(rse, assay, col = 1) unname(assays(rse)[[assay]][, col])
  test_plp_fxn <- function(...) {
    pileup_sites(bamfn, fafn, region = "SSR3:440-450",
               param = FilterParam(min_base_quality = 0,
                                         ...))}

  a <- test_plp_fxn()
  b <- test_plp_fxn(ftrim_5p = 0.10)
  d <- test_plp_fxn(ftrim_5p = 0.10, ftrim_3p = 0.10)
  e <- test_plp_fxn(ftrim_3p = 0.10)
  ex1 <- c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L)
  ex2 <- c(1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L)

  expect_equal(vec(a, "nRef") -  vec(b, "nRef"), ex1)
  expect_equal(vec(b, "nRef") - vec(d, "nRef"), ex2)
  expect_equal(vec(a, "nRef") - vec(e, "nRef"), ex2)

  f <- test_plp_fxn(ftrim_5p = 0.99, ftrim_3p = 0.99)
  expect_equal(length(f), 0)
  expect_error(test_plp_fxn(ftrim_3p = 1.10))
  expect_error(test_plp_fxn(ftrim_5p = -0.10))
})

fa <- scanFa(fafn)
hp_matches <- vmatchPattern(strrep("A", 6), fa) |> as("GRanges")

test_that("pileup check homopolymer filter", {
  a <- pileup_sites(bamfn, fafn, param = FilterParam(homopolymer_len = 6))
  b <- pileup_sites(bamfn, fafn)
  expect_false(identical(a, b))

  expect_equal(length(queryHits(findOverlaps(hp_matches, a))), 0)
  expect_equal(length(queryHits(findOverlaps(hp_matches, b))), 120)
})

test_that("filtering for splicing events works", {
  # get sites near splicing events
  splices <- GRanges(coverage(junctions(GenomicAlignments::readGAlignmentPairs(bamfn))))
  splices <- splices[splices$score > 0]
  rse <- pileup_sites(bamfn, fafn)
  rse <- rse[assay(rse, "ALT")[, 1] != "-"]

  # 3 sites, 1 is from all non-spliced reds, 2 sites are all from spliced reads,
  # (SPCS3  227, 347, 348)
  sites_near_splices <- rse[queryHits(findOverlaps(rse, splices, maxgap = 4))]

  rse <- pileup_sites(bamfn, fafn, param = FilterParam(splice_dist = 5))
  rse <- rse[assay(rse, "ALT")[, 1] != "-"]
  bsites_near_splices <- rse[queryHits(findOverlaps(rse,
    splices,
    maxgap = 4
  ))]
  expect_equal(nrow(sites_near_splices), 3)
  expect_equal(nrow(bsites_near_splices), 1)
})


test_that("filtering for indel events works", {
  # get sites near splicing events
  reads <- GenomicAlignments::readGAlignments(bamfn, use.names = T)
  cig_ops <- cigarRangesAlongReferenceSpace(cigar(reads), ops = "D", pos = start(reads))
  refs <- seqnames(reads)[elementNROWS(cig_ops) > 0]
  cig_ops <- cig_ops[elementNROWS(cig_ops) > 0]
  names(cig_ops) <- refs
  cig_pos <- unique(GRanges(cig_ops))

  rse <- pileup_sites(bamfn, fafn,
                    param = FilterParam(min_base_quality = 1,
                                              only_keep_variants = TRUE))

  # SSR3 387 variant has 1 read with variant with indel 4 bp away
  # SRR3 388 variant has 1 read with variant 3bp away from indel
  sites_near_indels <- rse[queryHits(findOverlaps(rse, cig_pos))]

  rse_b <- pileup_sites(bamfn, fafn, param = FilterParam(
    min_base_quality = 1,
    indel_dist = 5,
    only_keep_variants = TRUE
  ))
  bsites <- rse_b[queryHits(findOverlaps(rse_b, cig_pos))]

  expect_equal(nrow(sites_near_indels), 3)
  # lose SSR3 387 site after filtering
  expect_equal(nrow(bsites), 2)
  # lose 1 read from SSR3 388 site after filtering
  expect_equal( assay(sites_near_indels, "nAlt")[2] -
                  assay(bsites, "nAlt")[1], 1)
})

fout <- tempfile(fileext = ".fa")
test_that("writing reads with mismatches works", {
  rse <- pileup_sites(bamfn, fafn, reads = fout)
  seqs <- Rsamtools::scanFa(fout)
  expect_false(any(duplicated(names(seqs))))
})


seqs <- Rsamtools::scanFa(fout)
ids <- unlist(lapply(str_split(names(seqs), "_"), function(x) paste0(x[1], "_", x[2])))
writeLines(ids, fout)
test_that("excluding reads with mismatches works", {
  rse <- pileup_sites(bamfn, fafn, region = "DHFR:513-513")
  rse2 <- pileup_sites(bamfn, fafn,
    region = "DHFR:513-513",
    bad_reads = fout
  )
  expect_true(nrow(rse) == 1)
  expect_true(nrow(rse2) == 1)
  expect_true(assay(rse2, "nAlt") == 0)

  rse <- pileup_sites(bamfn, fafn)
  rse2 <- pileup_sites(bamfn, fafn, bad_reads = fout)
  expect_true(nrow(rse2) > 0)
  expect_true(colSums(assay(rse, "nAlt")) -
                colSums(assay(rse2, "nAlt")) > 0)
})

test_that("filtering for read-level mismatches works", {
  rse <- pileup_sites(bamfn, fafn,
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE
    ),
    region = "SSR3:244-247"
  )
  expect_equal(nrow(rse), 2)

  rse <- pileup_sites(bamfn, fafn,
    region = "SSR3:244-247",
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE,
      max_mismatch_type = c(1, 1)
    )
  )
  expect_equal(nrow(rse), 0)
  expect_s4_class(rse, "RangedSummarizedExperiment")

  rse <- pileup_sites(bamfn, fafn,
    region = "SSR3:244-247",
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE,
      max_mismatch_type = c(0, 10)
    )
  )
  expect_equal(nrow(rse), 2)

  rse <- pileup_sites(bamfn, fafn,
    region = "SSR3:244-247",
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE,
      max_mismatch_type = c(0, 1)
    )
  )
  expect_equal(nrow(rse), 0)
  expect_s4_class(rse, "RangedSummarizedExperiment")

  rse <- pileup_sites(bamfn, fafn,
    region = "SSR3:244-247",
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE,
      max_mismatch_type = c(8, 0)
    )
  )
  expect_equal(nrow(rse), 2)

  rse <- pileup_sites(bamfn, fafn,
    region = "SSR3:244-247",
    param = FilterParam(
      min_base_quality = 10,
      only_keep_variants = TRUE,
      max_mismatch_type = c(1, 0)
    )
  )
  expect_equal(nrow(rse), 0)
  expect_s4_class(rse, "RangedSummarizedExperiment")
})


test_that("parallel processing works", {
  rse <- pileup_sites(bamfn, fafn)
  rse_serial <- pileup_sites(bamfn, fafn, BPPARAM = SerialParam())
  expect_true(identical(rse_serial, rse))

  if (.Platform$OS.type != "windows") {
    rse_mc <- pileup_sites(bamfn, fafn, BPPARAM = MulticoreParam(workers = 2))
    expect_true(identical(rse_mc, rse))
  } else {
    rse_sp <- pileup_sites(bamfn, fafn, BPPARAM = SnowParam(workers = 2))
    expect_true(identical(rse_sp, rse))
  }
})


test_that("limiting chromosomes works", {
  chroms_to_query <- names(scanBamHeader(bamfn)[[1]]$targets)
  rse <- pileup_sites(bamfn, fafn)
  rse_2 <- pileup_sites(bamfn, fafn, chroms = chroms_to_query)
  expect_true(identical(rse_2, rse))

  rse_2 <- pileup_sites(bamfn, fafn, chroms = chroms_to_query[1])
  rse <- rse[seqnames(rse) == chroms_to_query[1]]
  expect_true(identical(rse_2, rse))

  if (.Platform$OS.type != "windows") {
    rse_mc <- pileup_sites(bamfn, fafn,
      chroms = chroms_to_query[1:2],
      BPPARAM = MulticoreParam(workers = 2)
    )
    expect_true(all(unique(seqnames(rse_mc)) == chroms_to_query[1:2]))
  }
})

unlink(c(bout, fout))


####### Other bam files

cbbam <- system.file("extdata", "5k_neuron_mouse_xf25_1pct_cbsort.bam", package = "raer")
tmp <- tempfile()
sort_cbbam <- sortBam(cbbam, tmp)
idx <- indexBam(sort_cbbam)

cbfa <- system.file("extdata", "mouse_tiny.fasta", package = "raer")

test_that("unsorted bam file fails", {
  expect_error(pileup_sites(cbbam, cbfa))
})

test_that("single end libary types are respected", {
  fp <- FilterParam(library_type = "genomic-unstranded")
  rse <- pileup_sites(sort_cbbam, cbfa, param = fp)
  expect_true(all(strand(rse) == "+"))
})

test_that("poor quality reads get excluded with read_bqual", {
  res_no_filter <- pileup_sites(bamfn, fafn,
    param = FilterParam(
      min_base_quality = 0L,
      library_type = c("genomic-unstranded")
    )
  )

  bqual_cutoff <- c(0.25, 20.0)
  bam_recs <- readGAlignments(bamfn,
    param = ScanBamParam(what = c("qname", "qual"))
  )
  strand(bam_recs) <- "+"
  bad_reads <- lapply(
    as(mcols(bam_recs)$qual, "IntegerList"),
    function(x) (sum(x < bqual_cutoff[2]) / length(x)) >= bqual_cutoff[1]
  )
  bad_reads <- bam_recs[unlist(bad_reads)]
  bad_reads <- subsetByOverlaps(bad_reads, res_no_filter)

  res <- pileup_sites(bamfn, fafn,
    param = FilterParam(
      read_bqual = bqual_cutoff,
      min_base_quality = 0L,
      library_type = c("genomic-unstranded")
    )
  )

  res_no_filter <- subsetByOverlaps(res_no_filter, bad_reads)
  res <- subsetByOverlaps(res, bad_reads)

  n_diff <- (assay(res_no_filter, "nRef") + assay(res_no_filter, "nAlt") -
               (assay(res, "nRef") + assay(res, "nAlt")))

  # we don't count deletions as coverage
  cov <- coverage(bad_reads, drop.D.ranges = TRUE)
  cov <- cov[cov != 0]
  cov <- as.integer(unlist(cov))

  expect_equal(as.vector(n_diff), cov)
})

test_that("sites within short splice overhangs can be excluded", {
  # 5 spliced reads, 2 unspliced
  fp <- FilterParam(
    library_type = "genomic-unstranded",
    min_splice_overhang = 0,
    min_base_quality = 0
  )
  res <- pileup_sites(bamfn, fafn, param = fp, region = "SPCS3:220-220")
  expect_equal(length(res), 1)
  expect_equal(as.vector(assay(res, "nRef")), 7)

  # excludes 1 read 5' of splice site (SRR5564269.2688337, cigar 51M120N100M)
  fp <- FilterParam(
    library_type = "genomic-unstranded",
    min_splice_overhang = 52,
    min_base_quality = 0
  )
  res <- pileup_sites(bamfn, fafn, param = fp, region = "SPCS3:220-220")
  expect_equal(length(res), 1)
  expect_equal(as.vector(assay(res, "nRef")), 6)

  fp <- FilterParam(
    library_type = "genomic-unstranded",
    min_splice_overhang = 0,
    min_base_quality = 0
  )
  res <- pileup_sites(bamfn, fafn, param = fp, region = "SPCS3:350-350")
  expect_equal(length(res), 1)
  expect_equal(as.vector(assay(res, "nRef")), 5)

  # excludes 1 read 3' of splice site (SRR5564269.224850, cigar 73M120N78M)
  fp <- FilterParam(
    library_type = "genomic-unstranded",
    min_splice_overhang = 79,
    min_base_quality = 0
  )
  res <- pileup_sites(bamfn, fafn, param = fp, region = "SPCS3:350-350")
  expect_equal(length(res), 1)
  expect_equal(as.vector(assay(res, "nRef")), 4)

  # excludes all 5 spliced reads
  fp <- FilterParam(
    library_type = "genomic-unstranded",
    min_splice_overhang = 101,
    min_base_quality = 0
  )
  res <- pileup_sites(bamfn, fafn, param = fp, region = "SPCS3:350-350")
  expect_equal(length(res), 0)
})
unlink(c(tmp, sort_cbbam, idx))

test_that("chroms missing from fasta file will not be processed", {
  mfafn <- raer_example("mouse_tiny.fasta")
  expect_error(expect_warning(pileup_sites(bamfn, mfafn, sites)))
})

test_that("excluding multiallelics works",{
  # DHFR_299_- has "A ->G,T"
  res_no_filter <- pileup_sites(bam2fn, fafn, region = "DHFR:299-299")
  expect_equal(assay(res_no_filter, "ALT")[1, 1], "T,G")

  fp <- FilterParam(report_multiallelic = FALSE)
  res <- pileup_sites(bam2fn, fafn, param = fp,  region = "DHFR:299-299")
  expect_equal(nrow(res), 0)

  fp <- FilterParam(report_multiallelic = FALSE)
  res <- pileup_sites(bam2fn, fafn, param = fp)
  n_og <- nrow(pileup_sites(bam2fn, fafn))
  expect_equal(nrow(res), n_og - 1)


  fp <- FilterParam(min_allelic_freq = 0.05)
  res <- pileup_sites(bam2fn, fafn, param = fp,  region = "DHFR:299-299")
  expect_equal(nrow(res), 1)
  expect_equal(assay(res, "ALT")[1, 1], "G")

  # ALT is first assay
  expect_equal(assays(res)[-1], assays(res_no_filter)[-1])

  res <- pileup_sites(bam2fn, fafn, param = fp)
  expect_equal(nrow(res), n_og)

  # the combination works,
  # e.g. excludes low frequency variant first, then checks if site is multiallelic
  fp <- FilterParam(min_allelic_freq = 0.05, report_multiallelic = FALSE)
  res <- pileup_sites(bam2fn, fafn, param = fp,  region = "DHFR:299-299")
  expect_equal(nrow(res), 1)
  expect_equal(assay(res, "ALT")[1, 1], "G")


  fp <- FilterParam(min_allelic_freq = 0.0001)
  res <- pileup_sites(bam2fn, fafn, param = fp,  region = "DHFR:299-299")
  expect_equal(nrow(res), 1)
  expect_equal(assay(res, "ALT")[1, 1], "T,G")

})






