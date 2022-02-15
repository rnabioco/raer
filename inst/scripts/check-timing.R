library(raer)

# ~8 Gb bam file
remote_bam <- "http://amc-sandbox.ucdenver.edu/User33/raer/SRR5564277_dedup_sorted.bam"
remote_bam <- "http://amc-sandbox.ucdenver.edu/User33/raer/SRR5564269_dedup_sorted.bam"
remote_fa <- "http://amc-sandbox.ucdenver.edu/User33/raer/Homo_sapiens.GRCh38.dna.primary_assembly.UCSC.fa"

# ~ 1 min for chr21
system.time(res <- get_pileup(bamfile = remote_bam, fafile = remote_fa, region = "chr21"))

# ~5 min for chr10
system.time(res <- get_pileup(bamfile = remote_bam, fafile = remote_fa, region = "chr10"))

# ~ 10 min for chr1
system.time(res <- get_pileup(bamfile = remote_bam, fafile = remote_fa, region = "chr1"))
