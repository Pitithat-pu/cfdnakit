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

This package provides functions for analyzing using shallow whole-genome sequencing data (~0.3X or more) of circulating cell-free DNA (cfDNA). The aims is to estimate circulating tumor DNA using its chracteristical short-fragmented cfDNA. The package extracts length of each cfDNA and assist the vistuallization of fragment-length distribution. A short-fragment ratio is calculated per non-overlapping fixed-sized bins. Genome-wide copy-number alteration estimated by the short-fragmented cfDNA . The ctDNA estimation score (CES) comprehensively estimate the circulating tumor DNA based on the short-fragment analysis.
}
\author{
Dr. rer. nat. Pitithat Puranachot
}


\keyword{ package cf-DNA}

\examples{
library(cfdnakit)
## Reading in a bamfile
sample_bamfile = system.file("extdata",
                             "ex.plasma.bam",
                package = "cfdnakit")
plasma_SampleBam = read_bamfile(sample_bamfile,
                         apply_blacklist = FALSE)

## Plot a fragment-length distribution of a sample
plot_fragment_dist(list("Plasma.Sample"=plasma_SampleBam))

## Plot a fragment-length distribution of two samples
control_RDS_file =
    system.file("extdata","BH01_CHR15.SampleBam.rds",
                package = "cfdnakit")
         ### Load example SampleBam of Healthy cfDNA
control_bins =
    readRDS(control_RDS_file)

comparing_list = list("Healthy.cfDNA"=control_bins,
                      "Patient.1"=plasma_SampleBam)
plot_fragment_dist(comparing_list)


## Derived and plot genome-wide short-fragment cfDNA
patient.SampleFragment =
  get_fragment_profile(plasma_SampleBam,
                       sample_id = "Patient.1")
plot_sl_ratio(patient.SampleFragment)


## Derived and plot normalized short-fragment cfDNA
PoN_rdsfile =  system.file(
                     "extdata",
                     "ex.PoN.rds",
                     package = "cfdnakit")
                     ## Loading example PoN data
PoN.profiles = readRDS(PoN_rdsfile)

sample_zscore =
  get_zscore_profile(patient.SampleFragment,
                     PoN.profiles)
sample_zscore_segment = segmentByPSCB(sample_zscore)
plot_transformed_sl(sample_zscore,sample_zscore_segment)

## Estimate circulating tumor DNA
calculate_CES_score(sample_zscore_segment)
}
