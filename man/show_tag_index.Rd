% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc-utils.R
\name{show_tag_index}
\alias{show_tag_index}
\title{Show tags stored in tag index}
\usage{
show_tag_index(bamfile)
}
\arguments{
\item{bamfile}{tag sorted bamfile, indexed with \code{\link[=build_tag_index]{build_tag_index()}}}
}
\value{
Character vector of tags
}
\description{
Show tags stored in tag index
}
\examples{
bam_fn <- raer_example("5k_neuron_mouse_xf25_1pct_cbsort.bam")
build_tag_index(bam_fn)
show_tag_index(bam_fn)

}
