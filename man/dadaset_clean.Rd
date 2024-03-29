% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dada2_processing_functions.R
\name{dadaset_clean}
\alias{dadaset_clean}
\title{Dadaset Clean}
\usage{
dadaset_clean(df, read_thresh = 100)
}
\arguments{
\item{df}{a merged dataframe i.e. taxonomy on rows and samples on columns.}

\item{read_thresh}{a numeric vector determining the number of reads per sample below which a sample is removed from the dataset. This cutoff is applied after filtering taxonomy. The default is 100.}
}
\value{
a merged dataframe which is free of unclassified taxonomy at the phylum level, and non-bacterial mitochondrial, and chloroplast sequences.
}
\description{
Removes taxonomy which is not classified as bacteria, or down to the phylum level, mitochondria, and chloroplast. Also removes samples which are below a designated read threshold after filtering.
}
\examples{
Theis_merged_clean <- dadaset_clean(df = Theis_merged, read_thresh = 100)
}
