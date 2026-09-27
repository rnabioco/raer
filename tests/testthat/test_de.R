pkgs <- c(
  "SummarizedExperiment",
  "SingleCellExperiment",
  "GenomicRanges",
  "Rsamtools"
)

msg <- lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")

rse_adar_ifn <- mock_rse()

test_that("calc_edit_frequency works", {
  # 6 sites had no coverage
  expect_message(rse <- calc_edit_frequency(rse_adar_ifn))

  expect_true("edit_freq" %in% assayNames(rse))
  expect_true("depth" %in% assayNames(rse))
  expect_true(all(assay(rse, "edit_freq") >= 0 & assay(rse, "edit_freq") <= 1))

  # depth already exists
  expect_message(calc_edit_frequency(rse))
  expect_error(calc_edit_frequency(
    rse,
    edit_from = "garbage-in",
    edit_to = "garbage-out"
  ))
})

test_that("make_de_object works", {
  rse_adar_ifn <- mock_rse()
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
  res <- find_de_sites(
    dse,
    condition_control = "WT",
    condition_treatment = "KO"
  )
  sig_sites <- rownames(res$sig_results)
  ed <- assay(rse, "edit_freq")[sig_sites, 1] -
    assay(rse, "edit_freq")[sig_sites, 4]
  expect_true(all(sign(ed) == sign(res$sig_results$logFC)))
})

test_that("find_scde_sites works", {
  sc_bamfn <- raer_example("5k_neuron_mouse_possort.bam")

  gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
  gr$REF <- c(rep("A", 4), "T")
  gr$ALT <- c(rep("G", 4), "C")

  cbs <- unique(scanBam(sc_bamfn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
  cbs <- na.omit(cbs)

  outdir <- tempdir()
  bai <- indexBam(sc_bamfn)

  fp <- FilterParam(library_type = "fr-second-strand")
  sce <- pileup_cells(sc_bamfn, gr, cbs, outdir, param = fp)

  set.seed(42)
  sce$clusters <- paste0("cluster_", sample(1:3, ncol(sce), replace = TRUE))
  res <- find_scde_sites(sce, "clusters")

  expect_true(is.list(res) || is(res, "List"))
  expect_setequal(
    names(res),
    unique(sce$clusters)
  )
  for (de_stats in as.list(res)) {
    expect_true(all(c("p.value", "dEF") %in% colnames(de_stats)))
    expect_true(all(de_stats$p.value >= 0 & de_stats$p.value <= 1, na.rm = TRUE))
  }
})
