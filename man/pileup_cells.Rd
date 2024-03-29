% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc-pileup.R
\name{pileup_cells}
\alias{pileup_cells}
\title{Generate base counts per cell}
\usage{
pileup_cells(
  bamfiles,
  sites,
  cell_barcodes,
  output_directory,
  chroms = NULL,
  umi_tag = "UB",
  cb_tag = "CB",
  param = FilterParam(),
  BPPARAM = SerialParam(),
  return_sce = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{bamfiles}{a path to a BAM file (for droplet scRNA-seq), or a vector of
paths to BAM files (Smart-seq2). Can be supplied as a character vector,
\link{BamFile}, or \link{BamFileList}.}

\item{sites}{a GRanges object containing sites to process. See examples for
valid formatting.}

\item{cell_barcodes}{A character vector of single cell barcodes to process.
If processing multiple BAM files (e.g. Smart-seq2), provide a character
vector of unique identifiers for each input BAM, to name each BAM file in the
output files.}

\item{output_directory}{Output directory for output matrix files. The
directory will be generated if it doesn't exist.}

\item{chroms}{A character vector of chromosomes to process. If supplied, only
sites present in the listed chromosomes will be processed}

\item{umi_tag}{tag in BAM containing the UMI sequence}

\item{cb_tag}{tag in BAM containing the cell-barcode sequence}

\item{param}{object of class \code{\link[=FilterParam]{FilterParam()}} which specify various filters to
apply to reads and sites during pileup. Note that the \code{min_depth} and
\code{min_variant_reads} parameters if set > 0 specify the number of reads
from any cell required in order to report a site. E.g. if
\code{min_variant_reads} is set to 2, then at least 2 reads (from any cell) must
have a variant in order to report the site. Setting \code{min_depth} and
\code{min_variant_reads} to 0 reports all sites present in the \code{sites} object.
The following options are not enabled for pileup_cells():
\code{max_mismatch_type}, \code{homopolymer_len}, and \code{min_allelic_freq}.}

\item{BPPARAM}{BiocParallel instance. Parallel computation occurs across
chromosomes.}

\item{return_sce}{if \code{TRUE}, data is returned as a SingleCellExperiment, if
\code{FALSE} a character vector of the output files, specified by
\code{outfile_prefix}, will be returned.}

\item{verbose}{Display messages}
}
\value{
Returns either a \link{SingleCellExperiment} or character vector of paths
to the sparseMatrix files produced. The \link{SingleCellExperiment} object is
populated with two assays, \code{nRef} and \code{nAlt}, which represent base counts
for the reference and alternate alleles. The \code{\link[=rowRanges]{rowRanges()}} will contain the
genomic interval for each site, along with \code{REF} and \code{ALT} columns. The
rownames will be populated with the format
\verb{site_[seqnames]_[position(1-based)]_[strand]_[allele]}, with \code{strand}
being encoded as 1 = +, 2 = -, and 3 = *, and allele being \code{REF} + \code{ALT}.

If \code{return_sce} is \code{FALSE} then a character vector of paths to the
sparseMatrix files (\code{barcodes.txt.gz}, \code{sites.txt.gz}, \code{counts.mtx.gz}),
will be returned. These files can be imported using \code{\link[=read_sparray]{read_sparray()}}.
}
\description{
This function processes scRNA-seq library to enumerate base
counts for Reference (unedited) or Alternate (edited) bases at specified
sites in single cells. \code{pileup_cells} can process droplet scRNA-seq
libraries, from a BAM file containing a cell-barcode and UMI, or well-based
libraries that do not contain cell-barcodes.

The \code{sites} parameter specifies sites to quantify. This must be a \link{GRanges}
object with 1 base intervals, a strand (+ or -), and supplemented with
metadata columns named \code{REF} and \code{ALT} containing the reference and
alternate base to query. See examples for the required format.

At each site, bases from overlapping reads will be examined, and counts of
each ref and alt base enumerated for each cell-barcode present. A single
base will be counted once for each UMI sequence present in each cell.
}
\examples{
library(Rsamtools)
library(GenomicRanges)
bam_fn <- raer_example("5k_neuron_mouse_possort.bam")

gr <- GRanges(c("2:579:-", "2:625:-", "2:645:-", "2:589:-", "2:601:-"))
gr$REF <- c(rep("A", 4), "T")
gr$ALT <- c(rep("G", 4), "C")

cbs <- unique(scanBam(bam_fn, param = ScanBamParam(tag = "CB"))[[1]]$tag$CB)
cbs <- na.omit(cbs)

outdir <- tempdir()
bai <- indexBam(bam_fn)

fp <- FilterParam(library_type = "fr-second-strand")
sce <- pileup_cells(bam_fn, gr, cbs, outdir, param = fp)
sce

# example of processing multiple Smart-seq2 style libraries

many_small_bams <- rep(bam_fn, 10)
bam_ids <- LETTERS[1:10]
 
# for unstranded libraries, sites and alleles should be provided on + strand 
gr <- GRanges(c("2:579:+", "2:625:+", "2:645:+", "2:589:+", "2:601:+"))
gr$REF <- c(rep("T", 4), "A")
gr$ALT <- c(rep("C", 4), "G")

fp <- FilterParam(
    library_type = "unstranded",
    remove_overlaps = TRUE
)

sce <- pileup_cells(many_small_bams,
    sites = gr,
    cell_barcodes = bam_ids,
    cb_tag = NULL,
    umi_tag = NULL,
    outdir,
    param = fp
)
sce

unlink(bai)

}
\seealso{
Other pileup: 
\code{\link{pileup_sites}()}
}
\concept{pileup}
