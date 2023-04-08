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

test_that("basic functionality works", {
    fp <- FilterParam(library_type = "fr-second-strand")
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

    expect_warning(pileup_cells(bam_fn, gr,
        c("non", "existent", "barcodes"),
        outdir,
        param = fp
    ))

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
    fp <- FilterParam(library_type = "fr-second-strand")
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_true(all(gr == rowRanges(sce)))

    # should drop 1 site with 0 reads
    fp <- FilterParam(
        library_type = "fr-second-strand",
        min_variant_reads = 1
    )
    sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
    expect_equal(nrow(sce), 4)
})


test_that("if multiple bams are supplied, treat each as a separate cell", {
    fp <- FilterParam(library_type = "fr-second-strand")
    sce <- pileup_cells(rep(bam_fn, 10),
        sites = gr,
        cell_barcodes = LETTERS[1:10],
        cb_tag = NULL,
        umi_tag = NULL,
        outdir, param = fp
    )
    expect_true(all(gr == rowRanges(sce)))
    expect_equal(ncol(sce), 10)
    expect_true(all(colnames(sce) == LETTERS[1:10]))
})

test_that("if multiple bams are supplied, require nbam # of barcodes", {
    fp <- FilterParam(library_type = "fr-second-strand")
    expect_error(pileup_cells(rep(bam_fn, 10),
        sites = gr,
        cell_barcodes = LETTERS[1:2],
        cb_tag = NULL,
        umi_tag = NULL,
        outdir, param = fp
    ))
})
