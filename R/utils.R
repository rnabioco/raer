#' Provide working directory for raer example files.
#'
#' @param path path to file
#'
#' @examples
#' raer_example("human.fasta")
#'
#' @return Character vector will path to an internal package file.
#' @export
raer_example <- function(path) {
    system.file("extdata", path, package = "raer", mustWork = TRUE)
}

# transpose a list
# https://stackoverflow.com/questions/30164803/fastest-way-to-transpose-a-list-in-r-rcpp
t_lst <- function(x) {
    split(unlist(x), sequence(lengths(x)))
}

check_tag <- function(x) {
    if (length(x) != 1 && nchar(x) != 2) {
        stop("supplied tag must by nchar of 2: ", x)
    }
}

chunk_vec <- function(x, n) {
    if (n == 1) {
        res <- list(`1` = x)
        return(res)
    }
    split(x, cut(seq_along(x), n, labels = FALSE))
}


# flatten top list, keeping names from inner list
unlist_w_names <- function(x) {
    nms <- unlist(lapply(x, names))
    res <- unlist(x, use.names = FALSE)
    names(res) <- nms
    res
}

#' Find regions with oligodT mispriming
#'
#' @description OligodT will prime at A-rich regions in an RNA. Reverse transcription
#' from these internal priming sites will install an oligodT sequence at the 3' end
#' of the cDNA. Sequence variants within these internal priming sites are enriched
#' for variants converting the genomic sequence to the A encoded by the oligodT primer.
#' Trimming poly(A) from the 3' ends of reads reduces but does not eliminate these signals
#'
#' This function will identify regions that are enriched for mispriming events. Reads
#' that were trimmed to remove poly(A) (encoded in the pa tag by 10x genomics) are
#' identified. The aligned 3' positions of these reads are counted, and sites passing
#' thresholds (at least 2 reads) are retained as possible sites of mispriming. Be default
#' regions 5 bases upstream and 20 bases downstream of these putative mispriming sites
#' are returned.
#'
#' @param bamfile path to bamfile
#' @param fafile path to fasta file
#' @param pos_5p distance 5' of mispriming site to define mispriming region
#' @param pos_3p distance 3' of mispriming site to define mispriming region
#' @param min_reads minimum required number of reads at a mispriming site
#' @param tag bam tag containing number of poly(A) bases trimmed
#' @param tag_values range of values required for read to be considered
#' @param n_reads_per_chunk number of reads to process in memory, see
#' [`Rsamtools::BamFile()`]
#' @param verbose if true report progress
#'
#' @returns A GenomicsRanges containing regions enriched for putative mispriming
#' events. The `n_reads` column specifies the number of polyA trimmed reads
#' overlapping the mispriming region. The `n_regions` column specifies the number
#' overlapping independent regions found in each chunk (dictated by `n_reads_per_chunk`).
#' The `A_freq` column indicates the frequency of A bases within the region.
#'
#' @importFrom Rsamtools ScanBamParam BamFile
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicRanges reduce resize trim
#' @importFrom IRanges grouplengths
#' @importFrom S4Vectors aggregate
#' @examples
#' bam_fn <- raer_example("5k_neuron_mouse_possort.bam")
#' fa_fn <- raer_example("mouse_tiny.fasta")
#' find_mispriming_sites(bam_fn, fa_fn)
#'
#' @export
find_mispriming_sites <- function(bamfile, fafile, pos_5p = 5, pos_3p = 20,
                                  min_reads = 2, tag = "pa", tag_values = 3:300,
                                  n_reads_per_chunk = 1e6, verbose = TRUE){

    if(pos_5p < 0 || pos_3p < 0) {
        cli::cli_abort("pos_5p and pos_3p must be positive integers")
    }

    tg_lst <- list(tag_values)
    names(tg_lst) <- tag
    sbp <- Rsamtools::ScanBamParam(tagFilter = tg_lst,tag = tag)
    bf <- Rsamtools::BamFile(bamfile, yieldSize = n_reads_per_chunk)
    open(bf)
    pa_pks <- GRanges()
    repeat {
        # should only return reads with pa tag set
        galn <- GenomicAlignments::readGAlignments(bf, param = sbp)

        if (length(galn) == 0) break
        gr <- as(galn, "GRanges")
        if(verbose) {
            s_ivl <- gr[1]
            e_ivl <- gr[length(gr)]
            message("working on ", s_ivl, " to ", e_ivl)
        }

        # count # of overlapping reads
        ans <- merge_pa_peaks(gr)
        pa_pks <- c(pa_pks, ans)
    }
    close(bf)
    # merge again, handle edge cases between yieldsizes
    ans <- GenomicRanges::reduce(pa_pks, with.revmap = TRUE)
    mcols(ans) <- S4Vectors::aggregate(pa_pks, mcols(ans)$revmap,
                                       mean_pal = mean(pa_pks$mean_pal),
                                       n_reads = sum(pa_pks$n_reads),
                                       drop = FALSE)

    # keep reads above threshold, slop, and merge adjacent misprimed regions
    ans <- ans[ans$n_reads >= min_reads]
    if(length(ans) == 0) {
        return(empty_mispriming_record())
    }

    ans <- GenomicRanges::resize(ans, pos_3p + width(ans))
    ans <- GenomicRanges::resize(ans, pos_5p + width(ans), fix = "end")
    ans <- GenomicRanges::trim(ans)

    res <- GenomicRanges::reduce(ans, with.revmap = TRUE)
    mcols(res) <- S4Vectors::aggregate(ans, mcols(res)$revmap,
                                       n_reads = sum(ans$n_reads),
                                       drop = FALSE)
    res$n_regions <- IRanges::grouplengths(res$grouping)
    res$grouping <- NULL
    res <- pa_seq_context(res, fafile)
    res
}

empty_mispriming_record <- function() {
    col_types <- list(
        n_reads = integer(),
        n_regions = integer(),
        A_freq = numeric()
    )
    df <- do.call(data.frame, col_types)
    gr <- GRanges(c(seqnames = NULL, ranges = NULL, strand = NULL))
    mcols(gr) <- df
    gr
}

merge_pa_peaks <- function(gr) {
    # get 3' end of read
    start(gr[strand(gr) == "+"]) <- end(gr[strand(gr) == "+"])
    end(gr[strand(gr) == "-"]) <- start(gr[strand(gr) == "-"])

    # merge and count reads within merged ivls
    ans <- GenomicRanges::reduce(gr, with.revmap = TRUE)
    mcols(ans) <- S4Vectors::aggregate(gr, mcols(ans)$revmap,
                                       mean_pal = mean(gr$pa),
                                       drop = FALSE)
    mcols(ans)$n_reads <- IRanges::grouplengths(ans$grouping)
    ans
}

#' @importFrom Rsamtools FaFile scanFa
#' @importFrom Biostrings letterFrequency reverseComplement
pa_seq_context <- function(gr, fafile){
    fa <- Rsamtools::FaFile(fafile)
    seqs <- Rsamtools::scanFa(fa, gr)
    seqs[strand(gr) == "-"] <- Biostrings::reverseComplement(seqs[strand(gr) == "-"])
    a_prop <- Biostrings::letterFrequency(seqs, "A") / width(gr)
    mcols(gr)$A_freq <- a_prop[, 1]
    gr
}

