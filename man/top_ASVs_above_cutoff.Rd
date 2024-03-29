% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dada2_processing_functions.R
\name{top_ASVs_above_cutoff}
\alias{top_ASVs_above_cutoff}
\title{Top ASVs Above Cutoff}
\usage{
top_ASVs_above_cutoff(
  ASV,
  TAX,
  META,
  cutoff,
  top,
  ASVs_on_Rows = TRUE,
  study,
  subject_col,
  subjects
)
}
\arguments{
\item{ASV}{numeric dataframe with ASVs on rows and samples on columns.}

\item{TAX}{taxonomy matching ASV/OTU dataset with taxa on rows.}

\item{META}{a metadata dataframe associated with samples on rows}

\item{cutoff}{percent relative abundance cutoff in decimal form above which ASVs will be listed.}

\item{top}{integer of ASVs above cutoff to return}

\item{ASVs_on_Rows}{boolean vector indicating if ASVs are on rows in \code{ASV}}

\item{study}{character vector name of study/title.}

\item{subject_col}{character vector of column name in \code{META} to group samples by.}

\item{subjects}{character vector of groups within \code{subject_col} to evaluate.}
}
\value{
a dataframe with ASV name, abundance in dataset, genus and species classifications, and exact sequence.
}
\description{
Determines ASVs above a certain cutoff and returns taxonomy and abundance data.
}
\examples{
Theis_top5 <- top_ASVs_above_cutoff(ASV = Theis_rare_ASV,
TAX = Theis_rare_TAX,
cutoff = 0.01,
top = 5,
ASVs_on_Rows = TRUE,
study = "Theis")
}
