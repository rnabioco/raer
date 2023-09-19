test_that("get_region works", {
    x <- get_region("chr1:1-20")
    expect_true(all(names(x) == c("chrom", "start", "end")))
    expect_true(x$chrom == "chr1" && x$start == 0 && x$end == 20)

    x <- get_region("chr1")
    expect_true(x$start == 0 && x$end == 2^31 - 1)
    expect_error(get_region("214124:!124:1:!24:chr1"))
})

test_that("filter mispriming works", {
    library(GenomicAlignments)
    bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
    fa_fn <- raer_example("mouse_tiny.fasta")

    expect_error(find_mispriming_sites("hello", fa_fn))
    expect_error(find_mispriming_sites(bam_fn, "world", verbose = FALSE))

    res <- find_mispriming_sites(bam_fn, fa_fn, min_reads = 1, verbose = FALSE)
    expect_true(length(res) == 1L)

    res <- find_mispriming_sites(bam_fn, fa_fn, min_reads = 100, verbose = FALSE)
    expect_true(length(res) == 0L)

    expect_error(find_mispriming_sites(bam_fn, fa_fn, tag = "nonsense"))
    expect_error(find_mispriming_sites(bam_fn, fa_fn, tag_values = "more-nonsense"))

    res <- find_mispriming_sites(bam_fn, fa_fn,
                                 min_reads = 1,
                                 pos_5p = 1,
                                 pos_3p = 1,
                                 verbose = FALSE)
    expect_true(width(res) == 3L)
})


test_that("correct_strand works", {
    bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
    fafn <- raer_example("human.fasta")
    fp <- FilterParam(library_type = "unstranded")
    rse <- pileup_sites(bamfn, fafn, param = fp)

    genes <- GRanges(c(
        "DHFR:200-400:+",
        "SPCS3:100-200:-"
    ))

    ans <- correct_strand(rse, genes)
    expect_true("-" %in% strand(ans))

    genic_rse <- subsetByOverlaps(rse, genes, ignore.strand = TRUE)
    expect_equal(nrow(ans), nrow(genic_rse))

    # sites overlapping ambiguous regions are discarded
    genes <- GRanges(c(
        "DHFR:200-400:+",
        "DHFR:200-400:-"
    ))
    ans <- correct_strand(rse, genes)
    expect_true(nrow(ans) == 0L)

    # - strand alts are complemented
    genes <- GRanges(c(
        "DHFR:200-400:-"
    ))

    ans <- correct_strand(rse, genes)
    alts <- assay(ans, "ALT")[, 1]
    alts <- alts[alts != "-" & as.vector(strand(ans)) ==  "-"]
    og_alts <- assay(rse[names(alts), ], "ALT")[, 1]
    expect_true(all(og_alts != alts))
    expect_true(all(nchar(og_alts) == nchar(alts)))

    # check all sites are retained if non-disjoint genic regions overlap all sites
    all_regions <- GRanges(seqinfo(rse))
    strand(all_regions) <- c("+", "-", "+")

    ans <- correct_strand(rse, all_regions)
    expect_equal(nrow(ans), nrow(rse))
    minus_rse <- subsetByOverlaps(ans, all_regions[2])
    expect_true(all(strand(minus_rse) == "-"))

    alts <- assay(minus_rse, "ALT")[, 1]
    alts <- alts[alts != "-" & as.vector(strand(minus_rse)) ==  "-"]
    og_alts <- assay(rse[names(alts), ], "ALT")[, 1]
    expect_true(all(comp_bases(alts) == og_alts))
})
