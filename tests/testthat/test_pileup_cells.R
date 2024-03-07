pkgs <- c("GenomicRanges", "Rsamtools", "SingleCellExperiment")
msg <- lapply(pkgs, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})

bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
indexBam(bam_fn)
gr <- GRanges(c(
    "2:579:-",
    "2:625:-",
    "2:645:-",
    "2:589:-",
    "2:601:-"
))
gr$REF <- c(rep("A", 4), "T")
gr$ALT <- c(rep("G", 4), "C")
gr <- sort(gr)

cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
cbs <- na.omit(cbs)

outdir <- tempdir()
on.exit(unlink(outdir))
fp <- FilterParam(library_type = "fr-second-strand")
sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)

test_that("basic functionality works", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_true(is(sce, "SingleCellExperiment"))
    expect_equal(dim(sce), c(5, 556))
    expect_true(all(colnames(sce) == cbs))
    expect_equal(nrow(sce), length(gr))
    fe <- all(c(
        "barcodes.txt.gz",
        "sites.txt.gz",
        "counts.mtx.gz"
    ) %in% dir(outdir))
    expect_true(fe)
    expect_true(is(assays(sce)$nRef, "dgCMatrix"))
    expect_true(sum(assays(sce)$nRef) > 0)
    expect_true(sum(assays(sce)$nAlt) > 0)

    # returns empty matrices
    sce <- pileup_cells(bam_fn, gr,
        c("non", "existent", "barcodes"),
        outdir,
        param = fp
    )
    expect_equal(sum(assay(sce, "nRef")) + sum(assay(sce, "nAlt")), 0L)

    expect_error(pileup_cells(bam_fn, as.data.frame(gr),
        cbs, outdir,
        fp = fp
    ))

    expect_error(pileup_cells(bam_fn, gr,
        1:10, outdir,
        param = fp
    ))
})

test_that("ranges are correct", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_true(all(gr == rowRanges(sce)))
})

test_that("depth filters work", {
    # should drop 1 site with 0 reads
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_variant_reads = 1
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 4)

    # should also drop 1 site with 0 reads
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_depth = 1
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 4)

    # drop 2 sites, 1 with only 7 variant reads and site with no depth
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_variant_reads = 10
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 3)

    # drops 2 sites, 1 with total depth of 230 and site with no depth
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_depth = 231
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 3)

    # also drop 1 site with 7 variant reads
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_depth = 231,
        min_variant_reads = 8
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 2)
})


test_that("if multiple bams are supplied, treat each as a separate cell", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    sce <- pileup_cells(rep(bam_fn, 10),
        sites = gr,
        cell_barcodes = LETTERS[1:10],
        cb_tag = NULL,
        umi_tag = NULL,
        outdir, param = fp
    )
    expect_true(all(gr == rowRanges(sce)))
    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), length(gr))
    expect_true(all(colnames(sce) == LETTERS[1:10]))
})

test_that("if multiple bams are supplied, make sure rows aren't duplicated", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 1)
    sce <- pileup_cells(rep(bam_fn, 10),
                        sites = gr,
                        cell_barcodes = LETTERS[1:10],
                        cb_tag = NULL,
                        umi_tag = NULL,
                        outdir, param = fp
    )
    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), 4)
    expect_true(all(colnames(sce) == LETTERS[1:10]))
})

test_that("if multiple bams are supplied, require nbam = # of barcodes", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    expect_error(pileup_cells(rep(bam_fn, 10),
        sites = gr,
        cell_barcodes = LETTERS[1:2],
        cb_tag = NULL,
        umi_tag = NULL,
        outdir, param = fp
    ))
})

test_that("output files are generated", {
    sce <- pileup_cells(bam_fn, gr, cbs, outdir)
    mtx_fns <- file.path(
        outdir,
        c(
            "counts.mtx.gz",
            "sites.txt.gz",
            "barcodes.txt.gz"
        )
    )
    expect_true(all(file.exists(mtx_fns)))
    are_gzipped <- lapply(mtx_fns, function(x) {
        fo <- file(x)
        xx <- summary(fo)$class == "gzfile"
        close(fo)
        xx
    }) |>
        unlist()
    expect_true(all(are_gzipped))
})

test_that("read_sparray reconstructs rse from output files", {
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = FilterParam(min_depth = 0))
    mtx_fns <- file.path(
        outdir,
        c(
            "counts.mtx.gz",
            "sites.txt.gz",
            "barcodes.txt.gz"
        )
    )
    sp_sce <- read_sparray(mtx_fns[1], mtx_fns[2], mtx_fns[3], "coordinate")
    expect_true(identical(sce, sp_sce))

    sp_sce <- read_sparray(mtx_fns[1], mtx_fns[2], mtx_fns[3], "index")
    expect_equal(rownames(sp_sce), as.character(1:5))
    expect_true(all(elementNROWS(rowRanges(sp_sce)) == 0))
})

test_that("BamFile and BamFileList input work", {
    bf <- BamFile(bam_fn)
    sce <- pileup_cells(bf, gr, cbs, outdir)
    expect_true(is(sce, "SingleCellExperiment"))
    expect_equal(dim(sce), c(4, 556))

    bfl <- BamFileList(bam_fn)
    sce <- pileup_cells(bfl, gr, cbs, outdir)
    expect_true(is(sce, "SingleCellExperiment"))
    expect_equal(dim(sce), c(4, 556))
})

test_that("custom indexes work", {
    bai_1 <- tempfile()
    bai_2 <- tempfile()
    file.copy(paste0(bam_fn, ".bai"), bai_1)
    file.copy(paste0(bam_fn, ".bai"), bai_2)

    bfl <- BamFileList(rep(bam_fn, 2), c(bai_2, bai_2))
    sce <- pileup_cells(bfl,
        sites = gr,
        cell_barcodes = LETTERS[1:2],
        cb_tag = NULL,
        umi_tag = NULL,
        outdir,
        param = FilterParam(min_depth = 0)
    )
    expect_true(all(gr == rowRanges(sce)))
    unlink(c(bai_1, bai_2))
})


# get larger set of sites to query
fa_fn <- raer_example("mouse_tiny.fasta")
bulkfp <- FilterParam(
    min_mapq = 255L,
    min_variant_reads = 1,
    min_allelic_freq = 0.01,
    only_keep_variants = TRUE,
    report_multiallelic = FALSE,
    library_type = "fr-second-strand",
    bam_flags = scanBamFlag(
        isSecondaryAlignment = FALSE,
        isSupplementaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE
    )
)
rse <- pileup_sites(bam_fn, fa_fn, chrom = "2", param = bulkfp)
rowData(rse)$ALT <- assay(rse, "ALT")[, 1]
sites <- rowRanges(rse)

test_that("rownames are not duplicated", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    sce <- pileup_cells(bam_fn, sites, cbs, outdir, param = fp)
    expect_false(any(duplicated(rownames(sce))))
})

test_that("multiple alleles per site can be queried", {
    fp <- FilterParam(library_type = "fr-second-strand", min_depth = 0)
    ol_sites <- GRanges(rep(c("2:261", "2:262", "2:279"), each = 2),
        strand = rep(c("+", "-"), 3),
        REF = c("G", "C", "C", "G", "G", "C"),
        ALT = c("T", "A", "T", "A", "T", "A")
    )
    sce <- pileup_cells(bam_fn, ol_sites, cbs, outdir, param = fp)
    expect_true(nrow(sce) == 6)
    expect_true(all(Matrix::rowSums(assay(sce, "nAlt")) > 0))
    expect_true(all(Matrix::rowSums(assay(sce, "nRef")) > 0))
})


mouse_sam_hdr <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:2\tLN:1115",
    "@SQ\tSN:6\tLN:400",
    "@SQ\tSN:11\tLN:493",
    "@SQ\tSN:8\tLN:396"
)

change_base <- function(aln, pos, new_base) {
    rec <- strsplit(aln, "\t")[[1]]
    seq <- rec[[10]]
    substr(seq, pos, pos) <- new_base
    rec[[10]] <- seq
    paste(rec, collapse = "\t")
}

change_qual <- function(aln, pos, new_qual, offset = 33) {
    rec <- strsplit(aln, "\t")[[1]]
    qual <- rec[[11]]
    substr(qual, pos, pos) <- intToUtf8(new_qual + offset)
    rec[[11]] <- qual
    paste(rec, collapse = "\t")
}

test_that("UMI consensus base selection works", {
    aln <- "A00836:50:HKF7VDSXX:4:1361:12337:2503\t16\t2\t247\t255\t14M315N77M\t*\t0\t0\tTTATTTATTTATTTTTTTGAGACAGGGTTTCTCTCTGTAGCCCTGGCTGTCCTGGAACTCACTCTGTAGACCAGGCTGGCCTCGAACTCAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGGAGTAGGCTATCT-1\tUB:Z:CTATGTTCGTGG\tRE:A:E\tNH:i:1\tHI:i:1\tnM:i:7\tCR:Z:AGGGAGTAGGCTATCT\tUR:Z:CTATGTTCGTGG\tAS:i:67\tCY:Z:FFFFFFFFFFFFFFFF\tUY:Z:FFFFFFFFFFFF\txf:i:0"
    ref_aln <- change_base(aln, 84, "A")

    tf <- tempfile(fileext = ".sam")
    writeLines(c(mouse_sam_hdr, ref_aln, aln), tf)
    tbam <- asBam(tf, overwrite = TRUE)
    tbai <- indexBam(tbam)
    on.exit(unlink(c(tf, tbam, tbai)))

    fp <- FilterParam(library_type = "fr-first-strand", min_depth = 0)
    ol_sites <- GRanges(c("2:20", "2:645", "2:646"),
        strand = "+",
        REF = "A",
        ALT = "G"
    )
    # 1 reads support ref 1 alt
    sce <- pileup_cells(tbam, ol_sites,
        "AGGGAGTAGGCTATCT-1",
        outdir,
        param = fp
    )
    expect_true(all(assay(sce, "nRef")[, 1] == c(0, 0, 1)))
    expect_true(all(assay(sce, "nAlt")[, 1] == c(0, 1, 0)))

    # 2 reads support ref 1 alt
    writeLines(c(mouse_sam_hdr, ref_aln, ref_aln, aln), tf)
    tbam <- asBam(tf, overwrite = TRUE)
    tbai <- indexBam(tbam)
    sce <- pileup_cells(tbam, ol_sites,
        "AGGGAGTAGGCTATCT-1",
        outdir,
        param = fp
    )
    expect_true(all(assay(sce, "nRef")[, 1] == c(0, 1, 1)))
    expect_true(all(assay(sce, "nAlt")[, 1] == c(0, 0, 0)))

    # 2 reads support ref, 1 support alt, alt has higher sum of base qualities
    ref_low_qual_aln <- change_qual(ref_aln, 84, 10)
    writeLines(c(mouse_sam_hdr, ref_low_qual_aln, ref_low_qual_aln, aln), tf)
    tbam <- asBam(tf, overwrite = TRUE)
    tbai <- indexBam(tbam)

    sce <- pileup_cells(tbam, ol_sites,
        "AGGGAGTAGGCTATCT-1",
        outdir,
        param = fp
    )
    expect_true(all(assay(sce, "nRef")[, 1] == c(0, 0, 1)))
    expect_true(all(assay(sce, "nAlt")[, 1] == c(0, 1, 0)))

    # no UMI tag disables consensus
    fp <- FilterParam(library_type = "fr-first-strand", min_base_quality = 10, min_depth = 0)
    sce <- pileup_cells(tbam, ol_sites,
        "AGGGAGTAGGCTATCT-1",
        umi_tag = NULL,
        outdir, param = fp
    )
    expect_true(all(assay(sce, "nRef")[, 1] == c(0, 2, 3)))
    expect_true(all(assay(sce, "nAlt")[, 1] == c(0, 1, 0)))

    # non ref/alt majority base not counted.
    other_aln <- change_base(aln, 84, "C")
    writeLines(c(mouse_sam_hdr, other_aln, other_aln, aln), tf)
    tbam <- asBam(tf, overwrite = TRUE)
    tbai <- indexBam(tbam)

    sce <- pileup_cells(tbam, ol_sites,
        "AGGGAGTAGGCTATCT-1",
        outdir,
        param = fp
    )
    expect_true(all(assay(sce, "nRef")[, 1] == c(0, 0, 1)))
    expect_true(all(assay(sce, "nAlt")[, 1] == c(0, 0, 0)))
})
