#' Download GSE99249 BAM files and related data
#'
#' @description  This function will download ~ 250 MB of data and store files in a
#' BiocFileCache.
#'
#' @param verbose report messages
#' @returns A named list with paths to BAM files, a FASTA file and a bed
#' file of known editing sites from hg38 chromosome 18.
#'
#' @family external-data
#'
#' @examples
#' \donttest{
#' download_GSE99249()
#' }
#'
#' @importFrom utils relist
#' @export
download_GSE99249 <- function(verbose = TRUE) {

    baseURL <- "https://raer-test-data.s3.us-west-2.amazonaws.com/GSE99249/"

    bam_fns <- c(
        "SRR5564260_dedup_sorted.bam",
        "SRR5564261_dedup_sorted.bam",
        "SRR5564269_dedup_sorted.bam",
        "SRR5564270_dedup_sorted.bam",
        "SRR5564271_dedup_sorted.bam",
        "SRR5564277_dedup_sorted.bam"
    )
    GSE99249_files <- list(
        bams = bam_fns,
        bai = paste0(bam_fns, ".bai"),
        fasta = "chr18.fasta.bgz",
        bed = "rediportal_hg38_chr18.bed.gz"
    )

    bfc <- .get_cache()
    fls <- unlist(GSE99249_files)
    fl_paths <- .add_files(fls, paste0(baseURL, fls), bfc, verbose)
    fl_paths <- utils::relist(fl_paths, GSE99249_files)
    fl_paths
}

#' Download NA12878 BAM files and related data
#'
#' @description This function will download ~100MB of data and store files in a
#' BiocFileCache.
#' @param verbose report messages
#'
#' @returns A named list with paths to an RNA-seq and WGS BAM file, and a FASTA file
#' from hg38 chromosome 4.
#'
#' @family external-data
#'
#' @examples
#' \donttest{
#' download_NA12878()
#' }
#'
#' @export
download_NA12878 <- function(verbose = TRUE) {

    baseURL <- "https://raer-test-data.s3.us-west-2.amazonaws.com/NA12878/"

    bam_fns <- c(
        "ERR262996_dedup_chr4_sub.bam",
        "SRR1258218_Aligned.sorted.dedup_sub.bam"
    )
    NA12878_files <- list(
        bams = bam_fns,
        bai = paste0(bam_fns, ".bai"),
        fasta = "hg38_chr4.fa.bgz",
        snps = "chr4snps.bed.gz",
        rmsk = "rmsk_hg38.tsv.gz"
    )

    bfc <- .get_cache()
    fls <- unlist(NA12878_files)
    fl_paths <- .add_files(fls, paste0(baseURL, fls), bfc, verbose)
    fl_paths <- utils::relist(fl_paths, NA12878_files)
    fl_paths
}

#' Download 10x PMBC bam file and related data
#'
#' @description This function will download < 1 GB of data and store files in a
#' BiocFileCache.
#' @param verbose report messages
#'
#' @returns A named list with paths to bam file, fasta file,
#' bed file of editing_sites, and an .rds file with a
#' SingleCellExperiment
#'
#' @family external-data
#'
#' @examples
#' \donttest{
#' download_human_pbmc()
#' }
#'
#' @export
download_human_pbmc <- function(verbose = TRUE) {
    baseURL <- "https://raer-test-data.s3.us-west-2.amazonaws.com/10x_human_pbmc/"

    pbmc_files <- list(
        bam = "10k_PBMC_3p_nextgem_Chromium_X_intron_possorted_chr16_rp.bam",
        bai = "10k_PBMC_3p_nextgem_Chromium_X_intron_possorted_chr16_rp.bam.bai",
        edit_sites = "rediportal_chr16.bed.gz",
        sce = "sce.rds"
    )

    bfc <- .get_cache()
    fls <- unlist(pbmc_files)
    fl_paths <- .add_files(fls, paste0(baseURL, fls), bfc, verbose)
    fl_paths <- utils::relist(fl_paths, pbmc_files)
    fl_paths
}

#' @importFrom BiocFileCache bfcquery bfcadd bfcneedsupdate bfcdownload bfcrpath
.add_files <- function(fls, urls, bfc, verbose = TRUE) {
  bam_suffixes <- c(".bai", ".bam")
  stopifnot((length(fls) == length(urls)) && (length(fls) > 0))

  fls_paths <- unlist(lapply(seq_along(fls), function(i) {
    fl <- fls[i]
    fl_url <- urls[i]
    is_bam <- any(vapply(bam_suffixes, endsWith, x = fl, logical(1)))
    rid <- BiocFileCache::bfcquery(bfc, fl, "rname")$rid
    if (!length(rid)) {
        if(verbose) {
            cli::cli_alert_info("Downloading file: {fl}")
        }
        rid <- names(BiocFileCache::bfcadd(bfc,
                                           fl,
                                           fpath = fl_url,
                                           fname = ifelse(is_bam,
                                                          "exact",
                                                          "unique")))
    }
    if (!isFALSE(any(BiocFileCache::bfcneedsupdate(bfc, rid)))){
      BiocFileCache::bfcdownload(bfc, rid)
    }
    BiocFileCache::bfcrpath(bfc, rids = rid[1])
  }))
  fls_paths
}

#' @importFrom BiocFileCache BiocFileCache
#' @importFrom tools R_user_dir
.get_cache <- function() {
    cache <- tools::R_user_dir("raer", which = "cache")
    BiocFileCache::BiocFileCache(cache)
}
