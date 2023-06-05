pkgs <- c(
    "GenomicRanges", "GenomicFeatures",
    "SummarizedExperiment", "rtracklayer"
)

msg <- lapply(pkgs, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})


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
        res <- annot_snps(se, SNPlocs.Hsapiens.dbSNP144.GRCh38)
        expect_true(is(res, "RangedSummarizedExperiment"))
        expect_true("RefSNP_id" %in% colnames(mcols(res)))

        if (require(BSgenome.Hsapiens.NCBI.GRCh38)) {
            res <- annot_snps(se,
                dbsnp = SNPlocs.Hsapiens.dbSNP144.GRCh38,
                genome = BSgenome.Hsapiens.NCBI.GRCh38
            )
            expect_true(is(res, "RangedSummarizedExperiment"))
            expect_true(all(c("RefSNP_id", "snp_ref_allele", "snp_alt_alleles") %in%
                colnames(mcols(res))))
        }
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


test_that("filter_multiallelic works", {
    data(rse_adar_ifn)
    x <- sum(grepl(",", assay(rse_adar_ifn, "ALT")))

    expect_message(rse <- filter_multiallelic(rse_adar_ifn))
    xx <- sum(grepl(",", assay(rse, "ALT")))

    expect_true(x > xx && xx == 0)
    expect_true("ALT" %in% names(rowData(rse)))
})

mock_txdb <- function() {
    gr <- GRanges(c(
        "DHFR:310-330:-",
        "DHFR:410-430:-",
        "SSR3:100-155:-",
        "SSR3:180-202:-"
    ))
    gr$source <- "raer"
    gr$type <- "exon"
    gr$source <- NA
    gr$phase <- NA_integer_
    gr$gene_id <- c(1, 1, 2, 2)
    gr$transcript_id <- rep(c("1.1", "2.1"), each = 2)
    makeTxDbFromGRanges(gr)
}

test_that("get_splice_sites and filtering works", {
    txdb <- mock_txdb()
    spl_sites <- get_splice_sites(txdb,
        slop = 3
    )
    expect_true(all(width(spl_sites) == 6))
    expect_error(get_splice_sites(spl_sites))

    data(rse_adar_ifn)
    expect_message(rse <- filter_splice_variants(rse_adar_ifn, txdb))
    n_removed <- nrow(subsetByOverlaps(rse_adar_ifn, rse, invert = TRUE))
    expect_equal(n_removed, 5L)
    expect_true("DHFR_328_-" %in% rownames(rse_adar_ifn))
    expect_false("DHFR_328_-" %in% rownames(rse))
})

test_that("removing clustered variants works", {
    txdb <- mock_txdb()

    data(rse_adar_ifn)
    expect_message(rse <- filter_multiallelic(rse_adar_ifn))

    gr <- rowRanges(rse)
    nd <- distanceToNearest(gr)
    gr1 <- gr[mcols(nd)$distance > 100]
    expect_message(gr2 <- filter_clustered_variants(rse, txdb))
    expect_true(identical(rowRanges(gr2), gr1))

    # check that transcript based removal works
    # DHFR_328 -> DHFR_423 overlaps 310-330 -> 410-430 splice
    ex <- exons(txdb)
    rse_ex <- subsetByOverlaps(rse, ex)

    expect_message(
        gr2 <- filter_clustered_variants(rse_ex,
            txdb,
            regions = "transcript",
            variant_dist = 20
        )
    )

    expt <- rowRanges(subsetByOverlaps(rse_ex, gr2, invert = T))
    ds <- distanceToNearest(mapToTranscripts(expt, txdb,
        extractor.fun = exonsBy
    ))
    expect_true(all(mcols(ds)$distance < 20))

    expect_message(
        gr2 <- filter_clustered_variants(rse_ex, txdb,
            variant_dist = 25
        )
    )
    expect_equal(length(gr2), 3L)
})

test_that("calc_confidence works", {
    data(rse_adar_ifn)
    rse <- calc_confidence(rse_adar_ifn)
    expect_true("confidence" %in% names(rowData(rse)))
    expect_true(identical(range(rowData(rse)$confidence), c(0, 1)))
})
