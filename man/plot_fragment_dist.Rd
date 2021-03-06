% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fragment_dist.R
\name{plot_fragment_dist}
\alias{plot_fragment_dist}
\title{Plot Fragment-length Distribution}
\usage{
plot_fragment_dist(readbam_list, maximum_length = 550, minimum_length = 20)
}
\arguments{
\item{readbam_list}{List; A list containing SampleBam object/objects from the read_bamfile function}

\item{maximum_length}{Int; Maximum length of fragment. cfDNA fragment longer than this value will not be considered; Default 550}

\item{minimum_length}{Int; Minimum length of fragment. cfDNA fragment shorter than this value will not be considered;  Default 20}
}
\value{
distribution plot
}
\description{
Plot Fragment-length Distribution
}
