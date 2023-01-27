#' Calculate the Adenosine Editing Index (AEI)
#'
#' @description The Adenosine Editing Index describes the magnitude of A-to-I
#'   editing in a sample. The index is a weighted average of editing events (G
#'   bases) observed at A positions. The vast majority A-to-I editing occurs in
#'   ALU elements in the human genome, and these regions have a high A-to-I
#'   editing signal compared to other regions such as coding exons. This
#'   function will perform pileup at specified repeat regions and return a
#'   summary AEI metric.
#'
#' @references Roth, S.H., Levanon, E.Y. & Eisenberg, E. Genome-wide
#' quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods
#' 16, 1131–1138 (2019). https://doi.org/10.1038/s41592-019-0610-9
#'
#' @param bam_fn bam file
#' @param fasta_fn fasta
#' @param alu_ranges GRanges or the name of a BEDfile with regions to query for
#'   calculating the AEI, typically ALU repeats. If a BED file is supplied it
#'   will not be filtered by the txdb option.
#' @param txdb A txdb object, if supplied, will be used to subset the alu_ranges
#'   to those found overlapping genes. Alternatively a GRanges object with gene
#'   coordinates.
#' @param snp_db either a SNPlocs package, GPos, or GRanges object. If supplied,
#'   will be used to exclude polymorphic positions prior to calculating the AEI.
#'   If `calc_AEI()` will be used many times, one could save some time by first
#'   identifying SNPs that overlap the supplied alu_ranges, and passing these as
#'   a GRanges to snp_db rather than supplying all known SNPs (see
#'   [get_overlapping_snps()]). Combined with using a bedfile for alu_ranges can
#'   also will save time.
#' @param filterParam object of class [FilterParam()] which specify various
#'   filters to apply to reads and sites during pileup.
#' @param BPPARAM A [BiocParallelParam] object for specifying parallel options
#'   for operating over chromosomes.
#' @param verbose report progress on each chromosome?
#'
#' @returns A named list with the AEI index computed for all allelic
#'   combinations. If correctly computed the signal from the A_G index should be
#'   higher than other alleles (T_C), which are most likely derived from noise
#'   or polymorphisms.
#'
#' @examples
#' suppressPackageStartupMessages(library(Rsamtools))
#' bamfn <- raer_example("SRR5564277_Aligned.sortedByCoord.out.md.bam")
#' fafn <- raer_example("human.fasta")
#' dummy_alu_ranges <- scanFaIndex(fafn)
#' calc_AEI(bamfn, fafn, dummy_alu_ranges)
#'
#' @importFrom BiocParallel bpstop bpmapply SerialParam
#' @importFrom GenomicFeatures genes
#' @importFrom rtracklayer export
#' @importFrom Rsamtools scanBamHeader
#' @importFrom IRanges subsetByOverlaps
#' @import S4Vectors
#' @import GenomicRanges
#'
#' @export
calc_AEI <- function(bam_fn,
                     fasta_fn,
                     alu_ranges = NULL,
                     txdb = NULL,
                     snp_db = NULL,
                     filterParam = FilterParam(),
                     BPPARAM = SerialParam(),
                     verbose = FALSE) {
  chroms <- names(Rsamtools::scanBamHeader(bam_fn)[[1]]$targets)

  if (length(bam_fn) != 1) {
    stop("calc_AEI only operates on 1 bam file at a time")
  }

  if (is.null(alu_ranges)) {
    warning(
      "querying the whole genome will be very ",
      "memory intensive and inaccurate.\n",
      "Consider supplying a GRanges object with ALU\n",
      "or related repeats for your species "
    )
  }

  genes_gr <- NULL
  tmp_files <- NULL
  alu_bed_fn <- NULL
  if (filterParam@library_type %in% c("unstranded", "genomic-unstranded")) {
    if (is.null(txdb)) {
      stop("txdb required for processing unstranded data")
    }
    filterParam@library_type == "genomic-unstranded"
    if (is(txdb, "TxDb")) {
      genes_gr <- suppressWarnings(GenomicFeatures::genes(txdb))
    } else {
      genes_gr <- txdb
    }
  }

  if (!is.null(alu_ranges)) {
    if (is(alu_ranges, "character")) {
      if (!file.exists(alu_ranges)) {
        stop("supplied alu ranges bedfile does not exist:\n",
          alu_ranges,
          call. = FALSE
        )
      }
      alu_bed_fn <- alu_ranges
    } else if (is(alu_ranges, "GRanges")) {
      alu_bed_fn <- tempfile(fileext = ".bed")
      tmp_files <- c(tmp_files, alu_bed_fn)
      if (!is.null(txdb)) {
        if (is(txdb, "TxDb")) {
          genes_gr <- suppressMessages(GenomicFeatures::genes(txdb))
        } else {
          genes_gr <- txdb
        }
        alu_ranges <- subsetByOverlaps(alu_ranges, genes_gr, ignore.strand = TRUE)
        alu_ranges <- reduce(alu_ranges)
      }
      chroms <- intersect(chroms, as.character(unique(seqnames(alu_ranges))))
      alu_ranges <- alu_ranges[seqnames(alu_ranges) %in% chroms]
      rtracklayer::export(alu_ranges, alu_bed_fn)
    } else {
      stop("unrecognized format for alu_ranges")
    }
  }

  snps <- NULL
  if (!is.null(snp_db)) {
    if (is(snp_db, "GRanges") || is(snp_db, "GPos")) {
      if (is(alu_ranges, "GRanges")) {
        snps <- subsetByOverlaps(snp_db, alu_ranges)
        snps <- split(snps, seqnames(snps))[chroms]
      } else {
        chroms <- intersect(chroms, as.character(unique(seqnames(snp_db))))
        snps <- split(snp_db, seqnames(snp_db))[chroms]
      }
    } else if (is(snp_db, "ODLT_SNPlocs")) {
      if (is(alu_ranges, "GRanges")) {
        snps <- snpsByOverlaps(snp_db, alu_ranges)
        snps <- split(snps, seqnames(snps))[chroms]
      } else {
        stop(
          "removing snps using a SNPloc package requires ",
          "alu_ranges to be supplied "
        )
      }
    } else {
      stop("unknown snpdb object type")
    }
  }
  if (is.null(snps)) {
    aei <- bpmapply(.calc_AEI_per_chrom,
      chroms,
      MoreArgs = list(
        bam_fn = bam_fn,
        fasta_fn = fasta_fn,
        alu_bed_fn = alu_bed_fn,
        filterParam = filterParam,
        snp_gr = NULL,
        genes_gr = genes_gr,
        verbose = verbose
      ),
      BPPARAM = BPPARAM,
      SIMPLIFY = FALSE
    )
  } else {
    if (length(chroms) != length(snps)) {
      stop("issue subsetting SNPdb and chromosomes")
    }
    aei <- bpmapply(.calc_AEI_per_chrom,
      chroms,
      snps,
      MoreArgs = list(
        bam_fn = bam_fn,
        fasta_fn = fasta_fn,
        alu_bed_fn = alu_bed_fn,
        filterParam = filterParam,
        genes_gr = genes_gr,
        verbose = verbose
      ),
      BPPARAM = BPPARAM,
      SIMPLIFY = FALSE
    )
  }
  bpstop(BPPARAM)

  names(aei) <- chroms
  aei_res <- lapply(seq_along(aei), function(x) {
    vals <- aei[[x]]
    id <- names(aei)[x]
    xx <- as.data.frame(t(do.call(data.frame, vals)))
    xx$allele <- rownames(xx)
    xx$chrom <- id
    rownames(xx) <- NULL
    xx
  })

  aei_res <- do.call(rbind, aei_res)
  aei_res <- split(aei_res, aei_res$allele)
  aei_res <- lapply(aei_res, function(x) 100 * (sum(x$alt) / (sum(x$ref) + sum(x$alt))))

  if (length(tmp_files) > 0) unlink(tmp_files)
  aei_res
}

.calc_AEI_per_chrom <- function(bam_fn,
                                fasta_fn,
                                alu_bed_fn,
                                chrom,
                                filterParam,
                                snp_gr,
                                genes_gr,
                                verbose) {
  if (verbose) {
    start <- Sys.time()
    message("\tworking on: ", chrom, " time: ", Sys.time())
  }
  filterParam@min_nucleotide_depth <- 1L
  filterParam@only_keep_variants <- FALSE

  plp <- pileup_sites(bam_fn,
    fafile = fasta_fn,
    bedfile = alu_bed_fn,
    chroms = chrom,
    filterParam = filterParam
  )

  if (!is.null(snp_gr) && !is.null(plp)) {
    plp <- subsetByOverlaps(plp, snp_gr,
      invert = TRUE,
      ignore.strand = TRUE
    )
  }

  if (filterParam@library_type == "genomic-unstranded") {
    plp <- correct_strand(plp, genes_gr)
  }

  if (verbose) {
    message("\tcompleted in : ", Sys.time() - start)
  }

  bases <- c("A", "T", "C", "G")
  var_list <- list()
  for (i in seq_along(bases)) {
    rb <- bases[i]
    other_b <- setdiff(bases, rb)
    j <- plp[rowData(plp)$Ref == rb]
    for (k in seq_along(other_b)) {
      ab <- other_b[k]
      id <- paste0(rb, "_", ab)
      n_alt <- sum(assays(j)[[paste0("n", ab)]][, 1])
      n_ref <- sum(assays(j)[[paste0("n", rb)]][, 1])
      var_list[[id]] <- c(
        alt = n_alt,
        ref = n_ref,
        prop = 0
      )
    }
  }
  var_list
}


#' Retrieve SNPs overlapping intervals
#'
#' @description This function will find SNPs overlapping supplied intervals
#'   using a SNPlocs package. The SNPs can be returned in memory (as GPos
#'   objects) or written to disk as a bed-file (optionally compressed).
#'
#' @param gr Intervals to query
#' @param snpDb A reference ot a SNPlocs database
#' @param output_file A path to an output file. If supplied the file can be
#'   optionally compressed by including a ".gz" suffix. If not supplied, SNPS
#'   will be returned as a [GenomicRanges::GPos] object
#'
#' @return GPos object containing SNPs overlapping supplied genomic intervals
#' @examples
#' if (require(SNPlocs.Hsapiens.dbSNP144.GRCh38)) {
#'   gr <- GRanges(rep("22", 10),
#'     IRanges(seq(10510077, 10610077, by = 1000)[1:10], width = 250),
#'     strand = "+"
#'   )
#'   get_overlapping_snps(gr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
#' }
#' @importFrom rtracklayer export
#' @importFrom BSgenome snpsByOverlaps
#'
#' @export
get_overlapping_snps <- function(gr,
                                 snpDb,
                                 output_file = NULL) {
  gr <- gr[seqnames(gr) %in% seqnames(snpDb)]

  # iterate through each contig, drop mcols (snpID) to reduce memory
  alu_snps <- vector("list", length = length(seqnames(snpDb)))
  for (i in seq_along(seqnames(snpDb))) {
    x <- seqnames(snpDb)[i]
    tmp_gr <- gr[seqnames(gr) == x]

    xx <- BSgenome::snpsByOverlaps(snpDb, tmp_gr)
    mcols(xx) <- NULL

    if (!is.null(output_file)) {
      rtracklayer::export(xx, output_file, append = TRUE, ignore.strand = TRUE)
    } else {
      alu_snps[[i]] <- xx
    }
  }
  if (!is.null(output_file)) {
    return(output_file)
  }
  unlist(as(alu_snps, "GRangesList"))
}

#' Apply strand correction using gene annotations
#'
#' @description Gene annotations are used to infer the likely strand of editing
#'   sites. This function will operate on unstranded datasets which have been
#'   processed using "genomic-unstranded" library type which reports variants
#'   with respect to the + strand for all sites. The strand of the editing site
#'   will be assigned the strand of overlapping features in the `genes_gr`
#'   object. Sites with no-overlap, or overlapping features with conflicting
#'   strands (+ and -) will be removed.
#'
#' @param rse RangedSummarizedExperiment object containing editing sites processed with
#'   "genomic-unstranded" setting
#' @param genes_gr GRanges object containing reference features to annotate the
#'   strand of the editing sites.
#'
#' @return RangedSummarizedExperiment object containing pileup assays,
#' with strand corrected based on supplied genomic intervals.
#'
#' @examples
#' suppressPackageStartupMessages(library("GenomicRanges"))
#'
#' bamfn <- raer_example("SRR5564269_Aligned.sortedByCoord.out.md.bam")
#' fafn <- raer_example("human.fasta")
#' fp <- FilterParam(library_type = "genomic-unstranded")
#' rse <- pileup_sites(bamfn, fafn, filterParam = fp)
#'
#' genes <- GRanges(c(
#'   "DHFR:200-400:+",
#'   "SPCS3:100-200:-",
#'   "SSR3:3-10:-",
#'   "SSR3:6-12:+"
#' ))
#'
#' correct_strand(rse, genes)
#'
#' @importFrom stringr str_count
#'
#' @export
correct_strand <- function(rse, genes_gr) {
  if (length(rse) == 0) {
    return(rse)
  }

  stopifnot(all(strand(rse) == "+"))
  #stopifnot(all(c("Ref", "Var", "nA", "nT", "nC", "nG") %in% names(mcols(gr))))

  genes_gr$gene_strand <- strand(genes_gr)
  rse <- annot_from_gr(rse, genes_gr, "gene_strand", ignore.strand = TRUE)

  # drop non-genic and multi-strand (overlapping annotations)
  rse <- rse[!is.na(rowData(rse)$gene_strand), ]
  rse <- rse[stringr::str_count(rowData(rse)$gene_strand, ",") == 0, ]

  flip_rows <- as.vector(strand(rse) != rowData(rse)$gene_strand)

  rowData(rse)$Ref[flip_rows] <- BASE_MAP[rowData(rse)$Ref[flip_rows]]

  flipped_variants <- vector(mode = "list", ncol(rse))
  to_flip <- assay(rse, "Var")[flip_rows, , drop = FALSE]
  flipped_variants <- apply(to_flip, c(1,2), function(x) {
    vapply(str_split(x, ","), function(y) paste0(unname(ALLELE_MAP[y]),
                                                 collapse = ","),
           character(1))
    })

  assay(rse, "Var")[flip_rows, ] <- flipped_variants


  # complement the nucleotide counts by reordering the assays
  assays_to_swap <- c("nA", "nT", "nC", "nG")
  og_order <- rownames(rse)
  sites_to_swap <- rse[flip_rows, , drop = FALSE]

  to_swap <- assays(sites_to_swap)[assays_to_swap]
  to_swap <- to_swap[c(2, 1, 4, 3)]
  names(to_swap) <- assays_to_swap
  assays(sites_to_swap)[assays_to_swap] <- to_swap

  rse <- rbind(rse[!flip_rows, ], sites_to_swap)
  rse <- rse[og_order, ]

  strand(rse) <- rowData(rse)$gene_strand
  rowData(rse)$gene_strand <- NULL
  rse
}


ALLELE_MAP <- c(
  TA = "CT",
  CA = "GT",
  GA = "AT",
  AT = "TC",
  CT = "GC",
  GT = "AC",
  AC = "TG",
  TC = "CG",
  GC = "AG",
  AG = "TA",
  TG = "CA",
  CG = "GA",
  `-` = "-"
)

BASE_MAP <- c(
  "A" = "T",
  "G" = "C",
  "C" = "G",
  "T" = "A",
  "N" = "N"
)
