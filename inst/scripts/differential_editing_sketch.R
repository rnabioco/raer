library(raer)
library(here)
library(Rsamtools)
library(BiocParallel)
library(data.table)
library(stringr)
library(SummarizedExperiment)
library(ggplot2)
library(edgeR)
library(ComplexHeatmap)

bams <- c(
  "SRR5564260",
  "SRR5564261",
  "SRR5564268",
  "SRR5564269",
  "SRR5564270",
  "SRR5564271",
  "SRR5564272",
  "SRR5564273",
  "SRR5564274",
  "SRR5564275",
  "SRR5564276",
  "SRR5564277"
)

res_dir <- here("results/GSE99249/")

# .rds file at /beevol/home/riemondy/dev/raerdata/results/GSE99249/merged_se.rds
if(!file.exists(file.path(res_dir, "merged_se.rds"))){

  # @ /beevol/home/riemondy/dev/raerdata/data/GSE99249/plps
  # scripts to generate @ /beevol/home/riemondy/dev/raerdata/results/GSE99249
  # 01_generate_known_plp.R + 01_generate_known_plp.sh

  plpdir <- here("data/GSE99249/plps")

  grs <- lapply(bams, function(x){
    res <- readRDS(file.path(plpdir, x, paste0(x, "_pileup.rds")))
    res <- suppressWarnings(unlist(as(res, "GRangesList")))
    sort(res)
  })
  names(grs) <- bams
  se <- create_se(grs)

  # replace missing values with 0, for now
  assay(se, "nA")[is.na(assay(se, "nA"))] <- 0
  assay(se, "nG")[is.na(assay(se, "nG"))] <- 0
  mcols(rowRanges(se))$Ref <- apply(assay(se, "Ref"),
                                    1,
                                    function(x) unique(x[!is.na(x)]))

  mcols(rowRanges(se))$edit_id <- paste0(seqnames(rowRanges(se)),
                                         "_",
                                         start(rowRanges(se)),
                                         "_",
                                         as.integer(strand(rowRanges(se))))
  rownames(se) <- mcols(rowRanges(se))$edit_id
  saveRDS(se, file.path(res_dir, "merged_se.rds"))
}

se <- readRDS(file.path(res_dir, "merged_se.rds"))
se_filtered <- se[mcols(rowRanges(se))$Ref == "A", ]

# compute editing frequencies (A -> G only)
assay(se_filtered, "edit_freq") <- assay(se_filtered, "nG") / (assay(se_filtered, "nA") + assay(se_filtered, "nG"))
assay(se_filtered, "edit_freq")[is.na(assay(se_filtered, "edit_freq"))] <- 0

# compute # of sites detected with 10 reads, > 0.01 edit freq
n_pass_filter <- colSums((assay(se_filtered, "edit_freq") > 0.01) &
                           ((assay(se_filtered, "nA") + assay(se_filtered, "nG")) >= 10))
colData(se_filtered)$n_sites <- n_pass_filter

edit_idx <- colSums(assay(se_filtered, "nG")) / (colSums(assay(se_filtered, "nA")) +
                                                   colSums(assay(se_filtered, "nG")))
colData(se_filtered)$edit_idx <- edit_idx

# add in metadata

# @ /beevol/home/riemondy/dev/raerdata/dbases/docs/PRJNA386593_download_log.txt
mdata <- fread(here("dbases/docs/PRJNA386593_download_log.txt"),
               data.table = FALSE)
mdata <- mdata[, c("run_accession", "library_name")]
mdata$library_name <-  str_replace(mdata$library_name, " biological replicate", "")
mdata$genotype <- str_extract(mdata$library_name , "^[A-Za-z0-9]+")
mdata$treatment <- str_extract(mdata$library_name , "Interferon beta|mock")
mdata$rep <- str_extract(mdata$library_name , "[0-9]$")
mdata$genotype_treatment <- str_c(mdata$genotype, " ", mdata$treatment)
mdata <- mdata[match(rownames(colData(se_filtered)), mdata$run_accession), ]
colData(se_filtered) <- cbind(colData(se_filtered), mdata)

cols <- RColorBrewer::brewer.pal(n = 9, "Set1")
p <- colData(se_filtered) %>%
  as.data.frame() %>%
  ggplot(aes(genotype_treatment, n_sites)) +
  geom_col(aes(fill = rep), position = position_dodge()) +
  scale_fill_manual(values = cols) +
  labs(x = NULL,
       y = "# of sites detected") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cowplot::save_plot(file.path(res_dir, "number_of_sites.pdf"), p, base_asp = 1)

p <- colData(se_filtered) %>%
  as.data.frame() %>%
  ggplot(aes(genotype_treatment, edit_idx)) +
  geom_col(aes(fill = rep), position = position_dodge()) +
  scale_fill_manual(values = cols) +
  labs(x = NULL,
       y = "Editing Index") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cowplot::save_plot(file.path(res_dir, "editing_index.pdf"), p, base_asp = 1)

# function to generate SummarizedExperiment for use with edgeR or DESeq2
# will generate a counts assay with a matrix formated with 2 columns per sample
prep_for_de <- function(se,
                        min_prop = 0.1,
                        max_prop = 0.9,
                        min_samples = 3){

  pass_cutoff <- (assay(se, "edit_freq") >= min_prop) & (assay(se, "edit_freq") <= max_prop)
  se <- se[rowSums(pass_cutoff) >= min_samples, ]

  ref <- assay(se, "nA")
  colnames(ref) <- paste0(colnames(ref), "_ref")
  alt <- assay(se, "nG")
  colnames(alt) <- paste0(colnames(alt), "_alt")
  res <- cbind(ref, alt)
  mdata <- colData(se)
  ref_mdata <- alt_mdata <- mdata
  rownames(ref_mdata) <- colnames(ref)
  ref_mdata$count <- "ref"
  rownames(alt_mdata) <- colnames(alt)
  alt_mdata$count <- "alt"
  mdata <- rbind(ref_mdata, alt_mdata)
  res <- SummarizedExperiment(assays = list(counts = res),
                              colData = mdata)
  res
}

deobj <- prep_for_de(se_filtered,
                     min_prop = 0.1,
                     max_prop = 0.9,
                     min_samples = 3)

# e.g. "SRR5564260" "SRR5564261" "SRR5564268"
sample <- deobj$sample

# e.g. "ADAR1KO" or "Wildtype"
condition <- deobj$genotype

# e.g. "ref" or "alt"
count <- deobj$count

design <- model.matrix(~0 + sample + condition:count)

# drop col to get full rank matrix
design <- design[,!grepl("countref", colnames(design))]


# set lib.sizes to 1 as doing comparisions within sample,
# perhaps not ideal?
dge <- DGEList(assay(deobj, "counts"), lib.size = rep(1, ncol(deobj)))
dge <- estimateDisp(dge)
fit <- glmFit(dge, design)

# need to fix column names in design
colnames(design) <- str_replace(colnames(design), ":", ".")
colnames(design) <- str_replace(colnames(design), " ", ".")

# test for difference between KO and WT
.cons <- "conditionADAR1KO.countalt-conditionWildtype.countalt"

# test for difference IFN-beta and mock
#.cons <- "conditionInterferon.beta.countalt-conditionmock.countalt"

.contrasts <- makeContrasts(
  contrasts = .cons,
  levels = design)

con <- glmLRT(fit, contrast = .contrasts)
con <- topTags(con, n = nrow(con))
con <- as.data.frame(con)
res <- con[order(con$FDR, con$PValue), ]
res <- res[res$FDR < 0.05, ]

to_plot <- rownames(res)[1:100]

edit_freqs <- assay(se_filtered, "edit_freq")[to_plot, ]

hmap <- Heatmap(edit_freqs,
                name = "Proportion edited",
                col = viridis::viridis(100),
                top_annotation = HeatmapAnnotation(genotype = se_filtered$genotype,
                                                   treatment = se_filtered$treatment,
                                                   col = list(
                                                     genotype =setNames(cols[1:2], unique(se_filtered$genotype)),
                                                     treatment = setNames(cols[3:4], unique(se_filtered$treatment)))),
                show_row_names = FALSE)

pdf(file.path(res_dir, "editing_heatmap.pdf"))
draw(hmap)
dev.off()
