library(Rsamtools)
bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
bam2fn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
dummy_alu_ranges <- scanFaIndex(fafn)

test_that("calc_aei works", {
  ko <- calc_AEI(bamfn, fafn, dummy_alu_ranges)
  wt <- calc_AEI(bam2fn, fafn, dummy_alu_ranges)
  expect_true(wt$A_G > ko$A_G)
  expect_true(names(wt)[which.max(unlist(wt))] == "A_G")
})


