pkgs <- c(
    "SummarizedExperiment"
)

msg <- lapply(pkgs, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")

data(rse_adar_ifn)

test_that("calc_edit_frequency works", {
    # 6 sites had no coverage
    expect_message(rse <- calc_edit_frequency(rse_adar_ifn))

    expect_true("edit_freq" %in% assayNames(rse))
    expect_true("depth" %in% assayNames(rse))
    expect_true(all(assay(rse, "edit_freq") >= 0 & assay(rse, "edit_freq") <= 1))

    # depth already exists
    expect_message(calc_edit_frequency(rse))
    expect_error(calc_edit_frequency(rse,
                                       edit_from = "garbage-in",
                                       edit_to = "garbage-out"))
})

test_that("make_de_object works", {
    data(rse_adar_ifn)
    expect_message(rse <- calc_edit_frequency(rse_adar_ifn))
    dse <- make_de_object(rse, min_samples = 1)
    expect_true("counts" %in% assayNames(dse))
    expect_equal(type(assay(dse)), "integer")
    expect_equal(2 * ncol(rse), ncol(dse))
    expect_equal(setdiff(c("ref", "alt"), dse$count), character(0))
})

test_that("find_de_sites works", {
    bams <- rep(c(bamfn, bam2fn), each = 3)
    sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
    names(bams) <- sample_ids

    fp <- FilterParam(only_keep_variants = TRUE)
    rse <- pileup_sites(bams, fafn, param = fp)
    rse$condition <- substr(rse$sample, 1, 2)

    expect_message(rse <- calc_edit_frequency(rse))
    dse <- make_de_object(rse)
    res <- find_de_sites(dse, condition_control = "WT", condition_treatment = "KO")
    sig_sites <- rownames(res$sig_results)
    ed <- assay(rse, "edit_freq")[sig_sites, 1] - assay(rse, "edit_freq")[sig_sites, 4]
    expect_true(all(sign(ed) == sign(res$sig_results$logFC)))
})
