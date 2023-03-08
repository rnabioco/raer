library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
bedfn <- raer_example("regions.bed")
sites <- import(bedfn)
res <- pileup_sites(bamfn, fafn)


test_that("annot_snps works", {
    if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
        gr <- GRanges(rep("22", 10),
            IRanges(
                seq(10510077,
                    10610077,
                    by = 1000
                )[1:10],
                width = 250
            ),
            strand = "+",
            seqinfo = seqinfo(SNPlocs.Hsapiens.dbSNP144.GRCh38)
        )
        gr <- keepStandardChromosomes(gr)
        genome(gr) <- "GRCh38.p2"
        res <- annot_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
        expect_true(is(res, "GRanges"))
        expect_true("RefSNP_id" %in% colnames(mcols(res)))

        gr2 <- gr
        expect_error(annot_snps(gr, gr2))

        se <- SummarizedExperiment(matrix(seq_along(gr)))
        rowRanges(se) <- gr
        rowData(se) <- gr
        res <- annot_snps(se, SNPlocs.Hsapiens.dbSNP144.GRCh38)
        expect_true(is(res, "RangedSummarizedExperiment"))
        expect_true("RefSNP_id" %in% colnames(mcols(res)))
    }
})


test_that("annot_from_gr works", {
    data(rse_adar_ifn)
    gr <- GRanges(rep(c("SSR3", "SPCS3"), c(5, 15)),
        IRanges(seq(1, 500, by = 25), width = 50),
        strand = "+"
    )

    gr$feature <- sample(1:100, size = 20)
    gr$id <- sample(LETTERS, size = 20)

    rse <- annot_from_gr(rse_adar_ifn, gr, c(feature_set = "feature", "id"))
    expect_true(all(c("feature_set", "id") %in% colnames(mcols(rse))))

    expect_error(annot_from_gr(rse_adar_ifn, gr, c(feature_set = "foo", "bar")))

    gr_res <- annot_from_gr(rowRanges(rse_adar_ifn), gr, c(feature_set = "feature", "id"))
    expect_true(all(c("feature_set", "id") %in% colnames(mcols(gr_res))))
})



test_that("calc_edit_frequency works", {
  data(rse_adar_ifn)
  rse <- calc_edit_frequency(rse_adar_ifn)
  expect_true("edit_freq" %in% assayNames(rse))
  expect_true("depth" %in% assayNames(rse))
})

test_that("prep_for_de works",{
  data(rse_adar_ifn)
  rse <- calc_edit_frequency(rse_adar_ifn)
  dse <- prep_for_de(rse, min_samples = 1)
  expect_true("counts" %in% assayNames(dse))
  expect_equal(type(assay(dse)), "integer")
  expect_equal(2 * ncol(rse), ncol(dse))
  expect_equal(setdiff(c("ref", "alt"), dse$count), character(0))
})

test_that("perform_de works",{
  bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
  bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
  fafn <- raer_example("human.fasta")

  bams <- rep(c(bamfn, bam2fn), each = 3)
  sample_ids <- paste0(rep(c("KO", "WT"), each = 3), 1:3)
  names(bams) <- sample_ids

  fp <- FilterParam(only_keep_variants = TRUE)
  rse <- pileup_sites(bams, fafn, param = fp)
  rse$condition <- substr(rse$sample, 1, 2)

  rse <- calc_edit_frequency(rse)
  dse <- prep_for_de(rse)
  res <- perform_de(dse, condition_control = "WT", condition_treatment = "KO")
  sig_sites <- rownames(res$sig_results)[1:5]
  ed <- assay(rse, "edit_freq")[sig_sites, 1] - assay(rse, "edit_freq")[sig_sites, 4]
  expect_true(all(sign(ed) == sign(res$sig_results$logFC[1:5])))
})
