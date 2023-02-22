library(raer)

ko_bam <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
wt_bam <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
fafn <- raer_example("human.fasta")
bedfn <- raer_example("regions.bed")

fp <- FilterParam(
  only_keep_variants = TRUE,
  library_type = "fr-first-strand",
  min_depth = 2
)

bams <- c(wt_bam, ko_bam)
names(bams) <- c('wt', 'adar1_ko')

rse_adar_ifn <- pileup_sites(
  bams,
  fafn,
  param = fp
)

usethis::use_data(rse_adar_ifn, overwrite = TRUE, compress = 'xz')
