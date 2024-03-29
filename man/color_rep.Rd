% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dada2_processing_functions.R
\name{color_rep}
\alias{color_rep}
\title{Color Replication}
\usage{
color_rep(colors, full_set, highlight, highlight_bool, nonhgrey)
}
\arguments{
\item{colors}{color palette as character vector of colors to be used.}

\item{full_set}{factor vector indicating sample group for each sample.}

\item{highlight}{character vector of all samples to be colored normally}

\item{highlight_bool}{boolean vector indicating which samples are in \code{highlight}.}

\item{nonhgrey}{boolean vector indicating whether non-highlighted samples should be grey (\code{TRUE}) or left blank (\code{FALSE})}
}
\value{
a new character vector containing the new colors by sample.
}
\description{
Auxillary function to \code{beta_div()} which allows for modification of plotted colors. In particular, it allows for certain points to be highlighted while others are colored grey or left blank.
}
