pkgs <- c("Rsamtools", "data.table", "Biostrings")
msg <- lapply(pkgs, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})


bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
bams <- c(bamfn, bam2fn)
names(bams) <- c("ko", "wt")
fafn <- raer_example("human.fasta")
mock_alu_ranges <- scanFaIndex(fafn)
set.seed(42)

# snps at all As in SPCS3
mock_snps <- GRanges("SPCS3", IRanges(matchPattern("A", scanFa(fafn)[["SPCS3"]])))
mock_snps_gpos <- as(mock_snps, "GPos")

mock_genes <- mock_alu_ranges
strand(mock_genes) <- c("-","+", "-")

test_that("calc_aei basic options work", {
    fp <- FilterParam(library_type = "fr-first-strand")
    aei <- calc_AEI(bams, fafn, mock_alu_ranges, param = fp)
    expect_equal(length(aei), 2L)
    expect_equal(dim(aei$AEI), c(2L, 12L))
    expect_true(all(rownames(aei$AEI) == names(bams)))
    ag_aei <- aei$AEI[, "A_G"]
    expect_true(ag_aei["wt"] > ag_aei["ko"])

    x <- aei$AEI_per_chrom
    setDT(x)
    pc_aei <- x[, 100 * (sum(alt) / sum(alt + ref)), by = .(allele, bam_file)]
    pc_aei <- dcast(pc_aei, bam_file ~ allele, value.var = "V1")
    pc_aei <- as.data.frame(pc_aei)
    rownames(pc_aei) <- pc_aei$bam_file
    pc_aei$bam_file <- NULL
    expect_true(identical(pc_aei, as.data.frame(aei$AEI)))

    aei <- calc_AEI(unname(bams), fafn, mock_alu_ranges)
    expect_true(all(rownames(aei$AEI) == unname(bams)))

    # drops all As on SSR3 chrom
    aei_snp <- calc_AEI(bams, fafn, mock_alu_ranges, snp_db = mock_snps, param = fp)
    pc <- aei_snp$AEI_per_chrom
    ssr3_aei <- pc[pc$chrom == "SPCS3" & startsWith(pc$allele, "A_"), ]
    expect_true(all(ssr3_aei[c("alt", "ref")] == 0))

    pc <- aei$AEI_per_chrom
    ssr3_aei <- pc[pc$chrom == "SPCS3" & startsWith(pc$allele, "A_"), ]
    expect_true(all(ssr3_aei$ref > 0))

    # GPos works
    aei_snp <- calc_AEI(bams, fafn, mock_alu_ranges, snp_db = mock_snps_gpos, param = fp)
    pc <- aei_snp$AEI_per_chrom
    ssr3_aei <- pc[pc$chrom == "SPCS3" & startsWith(pc$allele, "A_"), ]
    expect_true(all(ssr3_aei[c("alt", "ref")] == 0))

    fp <- FilterParam(library_type = "genomic-unstranded")
    expect_error(calc_AEI(bams, fafn, mock_alu_ranges, param = fp))

    us_aei <- calc_AEI(bams, fafn, mock_alu_ranges, mock_genes, param = fp)
    expect_true(all(us_aei$AEI[, "A_G"] > us_aei$AEI[, "T_C"]))

})
