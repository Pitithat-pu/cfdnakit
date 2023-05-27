\name{cfdnakit-package}
\alias{cfdnakit-package}
\alias{cfdnakit}
\docType{package}
\title{
\packageTitle{cfdnakit}
}
\description{
\packageDescription{cfdnakit}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{cfdnakit}
\packageIndices{cfdnakit}
This package provides functions for analyzing using shallow whole-genome sequencing data (~0.3X or more) of circulating cell-free DNA (cfDNA). The aims is to estimate circulating tumor DNA using its chracteristical short-fragmented cfDNA. The package extracts length of each cfDNA and assist the vistuallization of fragment-length distribution. A short-fragment ratio is calculated per non-overlapping fixed-sized bins. Genome-wide copy-number alteration estimated by the short-fragmented cfDNA . The ctDNA estimation score (CES) comprehensively estimate the circulating tumor DNA based on the short-fragment analysis.
}
\author{
\packageAuthor{cfdnakit}
Dr. rer. nat. Pitithat Puranachot
Lecturer at Chulabhorn Royal Academy, Bangkok, Thailand
Maintainer: \packageMaintainer{cfdnakit}
Dr. rer. nat. Pitithat Puranachot
Lecturer at Chulabhorn Royal Academy, Bangkok, Thailand
}


\keyword{ package cf-DNA}

\examples{
library(cfdnakit)
### Reading in a bamfile
plasma_SampleBam = read_bamfile("Path.to/bamfile",apply_blacklist = FALSE)
### Plot a fragment-length distribution of a sample
plot_fragment_dist(list("Plasma.Sample"=plasma_SampleBam))

### Plot a fragment-length distribution of two samples
comparing_list = list("Healthy.cfDNA"=control_bins,
                      "Patient.1"=sample_bins)
plot_fragment_dist(comparing_list)


### Derived and plot genome-wide short-fragment cfDNA
sample_profile =
  get_fragment_profile(sample_bins,
                       sample_id = "Patient.1")
plot_sl_ratio(patient.SampleFragment)


### Derived and plot normalized short-fragment cfDNA
sample_zscore =
  get_zscore_profile(patient.SampleFragment,
                     PoN.profiles)
sample_zscore_segment = segmentByPSCB(sample_zscore)
plot_transformed_sl(sample_zscore,sample_zscore_segment)


### Estimate circulating tumor DNA
calculate_CES_score(sample_zscore_segment)
}
