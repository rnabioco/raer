#format reditools output
library(data.table)
dat <- fread("~/rbi/src/raerdata/TABLE1_hg38.txt.gz", data.table = FALSE)
bed <- dat[c("Region", "Position", "Strand")]
colnames(bed) <- c("chrom", "start", "strand")
bed$end <- bed$start
bed$start <- bed$start - 1L
chroms_to_keep <- paste0("chr", 1:22)
bed <- bed[bed$chrom %in% chroms_to_keep, ]
bed <- bed[order(bed$chrom, bed$start, bed$end), ]
fwrite(bed[c("chrom", "start", "end", "strand")],
       "rediportal_hg38.bed",
       quote = F, col.names = F, sep = "\t")
