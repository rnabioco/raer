library(GenomicRanges)
library(Rsamtools)

bam_fn <- raer_example("5k_neuron_mouse_possort.bam")

gr <- GRanges(c("2:579:-",
                "2:625:-",
                "2:645:-",
                "2:589:-",
                "2:601:-"))
gr$ref <- c(rep("A", 4), "T")
gr$alt <- c(rep("G", 4), "C")
gr <- sort(gr)

cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
cbs <- na.omit(cbs)

outdir <- tempdir()
on.exit(unlink(outdir))

bai <- Rsamtools::indexBam(bam_fn)
on.exit(unlink(bai))

test_that("basic functionality works", {
  fp <- FilterParam(library_type = "fr-second-strand")
  sce <- pileup_cells(bam_fn, gr, cbs, outdir, fp = fp)
  expect_true(is(sce, "SingleCellExperiment"))
  expect_equal(dim(sce), c(5, 556))
  expect_true( all(colnames(sce) == cbs))
  expect_equal(nrow(sce), length(gr))
  fe <- all(c("barcodes.txt.gz",
              "sites.txt.gz",
              "counts.mtx.gz") %in% dir(outdir))
  expect_true(fe)
  expect_true(is(assays(sce)$nRef, "dgCMatrix"))
  expect_true(sum(assays(sce)$nRef) > 0)
  expect_true(sum(assays(sce)$nVar) > 0)

  expect_warning(pileup_cells(bam_fn, gr,
                              c("non", "existent", "barcodes"),
                              outdir, fp = fp))

  expect_error(pileup_cells(bam_fn, as.data.frame(gr),
                            cbs, outdir, fp = fp))
  expect_error(pileup_cells(bam_fn, gr,
                            1:10, outdir, fp = fp))
})

test_that("ranges are correct", {
  fp <- FilterParam(library_type = "fr-second-strand")
  sce <- pileup_cells(bam_fn, gr, cbs, outdir, fp = fp)
  expect_true(all(gr == rowRanges(sce)))

  # should drop 1 site with 0 reads
  fp <- FilterParam(library_type = "fr-second-strand",
                    min_variant_reads = 1)
  sce <- pileup_cells(bam_fn, gr, cbs, outdir, fp = fp)
  expect_equal(nrow(sce), 4)
})
