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

