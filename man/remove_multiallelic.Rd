% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_se.R
\name{remove_multiallelic}
\alias{remove_multiallelic}
\title{Remove multiallelic sites}
\usage{
remove_multiallelic(se)
}
\arguments{
\item{se}{SummarizedExperiment}
}
\description{
Sites with multiple variant bases will be removed from the
SummarizedExperiment. The rowData() will be updated to include a new column,
"Var", containing the variant allele detected at each site.
}
\examples{
example(create_se, echo = FALSE)
remove_multiallelic(se)

}
