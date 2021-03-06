% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dada2_processing_functions.R
\name{list_types}
\alias{list_types}
\title{Create Sample Type List}
\usage{
list_types(df, types_on_col = FALSE, split_by = "_")
}
\arguments{
\item{df}{dataset with sample names on columns or rows.}

\item{types_on_col}{if sample names are on columns this should be true, if on rows this should be false}

\item{split_by}{character to separate sample type name from unique identifier. Default is underscore.}
}
\value{
Returns list of sample types and their counts based on the first letter of every sample.
}
\description{
Crude way to create Sample Type Lists for later PCoA plotting, heatmapping, or barplotting analyses.
}
\examples{
Sample_Types <- list_types(merged_clean, types_on_col = TRUE, split_by = "_")
}
