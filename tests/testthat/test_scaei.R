pkgs <- c("GenomicRanges", "Rsamtools", "SingleCellExperiment")
msg <- lapply(pkgs, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})

bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
indexBam(bam_fn)
fa_fn <- raer_example("mouse_tiny.fasta")

# cell barcodes to query
cbs <- c("TGTTTGTTCCATCCGT-1", "CAACCAACATAATCGC-1", "TGGAACTCAAGCTGTT-1")

genes_gr <- GRanges(c(
    "2:100-400:-",
    "2:500-605:-",
    "2:600-680:+"
))

# alu intervals
alus_gr <-  GRanges(c(
    "2:110-380",
    "2:510-600",
    "2:610-670"
))

test_that("calc_scaei basic functions work", {
    sites <- get_scAEI_sites(fa_fn, genes_gr, alus_gr)

    s <- scanFa(fa_fn, sites)
    s[strand(sites) == "-"] <- reverseComplement(s[strand(sites) == "-"])
    expect_true(all(s == "A"))

    fp <- FilterParam(library_type = "fr-second-strand",
                      min_mapq = 255)
    res <- calc_scAEI(bam_fn, sites, cbs, fp)
    expect_equal(nrow(res), 3L)
    expect_true(is(res, "DataFrame"))

    res <- calc_scAEI(bam_fn, sites, cbs, fp, return_sce = TRUE)
    expect_true(is(res, "SingleCellExperiment"))
})

