# cfdnakit : an R package for fragmentation analysis of cfDNA and copy-number alteration calling

## Introduction

This package provides mostly basic functions to analysis cell-free DNA (cfDNA) with next-generation sequencing. The package focuses on extracting length of cfDNA, and genome-wide copy-number alteration estimated by the short-fragmented cfDNA using shallow whole-genome sequencing data (\~0.3X or more).

### Overview

The figure below shows the overview of the CNV-calling procedure. The amount of short-fragmented cfDNA per non-overlapping genomics windows are normalised and compared to a Panel-of-Normal (control). Segmentation is performed using the PSCBS package. A modified CPA score (adapted from [Raman, Lennart, et al. 2020](https://doi.org/10.1186/s13073-020-00735-4)) were calculated to estimate the copy number tumor burden.

<img src="https://raw.githubusercontent.com/wiki/Pitithat-pu/cfdnakit/images/cfdnakit_workflow.png" title="cfdnakit cnv calling workflow" alt="cfdnakit_workflow" width="600"/>

## Available functions

The scope of this R package is to analyse the length of cfDNA fragments. The package simplifies the process of extracting length of fragments from a BAM file and provides basic functions to explore this characteristic of cfDNA with low-coverage whole-genome sequencing data. Moreover, this package utilizies the quantity of short-fragmented cfDNA to infer copy-number alterations and estimate the percentage of tumor-derived cfDNA.

### Excluding ENCODE blacklist loci

It is recommended when analyzing genomic data to exclude sequencing reads locate within the ENCODE blacklist loci to assure the quality of the result. When using the GRCh37 as the reference in cfdnakit, a set of genomic regions including the ENCODE blacklist and centromere loci, provided by UCSC Genome Browser, were used. Users can introduce customized blacklist regions by creating a bed file or a tab-separated file where the first three columns are chromosome, start, and end position respectively. The future cfdnakit would be able to support blacklists of other reference genomes such as GRCh38 or GRCm38.

### Fragment length distribution and comparison

Package provides a single function to extract fragment length of cfDNA in the sample. Making a fragment-length distribution plot of multiple samples is easy. cfdnakit also extracted the short-fragment ratio representing the amount of short-fragmented cfdNA in the sample. It can be used for comparison between groups of sample (e.g. healthy vs patient) or for quality control inspection.

This plot shows the fragment-length distribution of cfDNA from a healthy individual cfDNA (red) and a patient-derived cfDNA (blue). The top-right legend report the leak length per individual sample as shown.

<img src="https://github.com/Pitithat-pu/cfdnakit/wiki/images/wiki/fragment_length_dis.png" title="cfdnakit fragment length distribution" alt="fragment_length_distribution_cfdnakit" width="550" height="420"/>

### GC and Mappability Bias Correction

LOESS regression model is created from the relation between the fragment count and the percent of GC per bin. The raw count per bin is deduced with the read count predicted by the model. Then, the values are added with the median of raw counts to bring back the range of values similar to the raw count. After correction for GC bias, the GC-corrected read counts are then corrected for mappability bias using a similar process. The mappability bias indicates the mapping capability of a genomic region to be mapped uniquely by sequencing reads.

### Sample Short-fragmented Ratio

An enrichment of short-fragmented cfDNA (\<150 base) were commonly found in tumor-derived cfDNA. This package extract a preliminary estimation of short-fragmented cfDNA by calculating Short-fragmented Ratio (S.L.Ratio) inferring to the proportion of short-fragmented cfDNA (100 to 150 by default) over long-fragmented cfDNA (151 to 250 by default).

![](https://latex.codecogs.com/svg.image?S.L.Ratio_%7Bsample%7D=%5Cfrac%7BN_%7BF.short%7D%7D%7BN_%7BF.long%7D%7D)

where ![](https://latex.codecogs.com/svg.image?N_%7BF.short%7D) is number of short fragments; ![](https://latex.codecogs.com/svg.image?N_%7BF.long%7D) is number of long fragments;

![](https://latex.codecogs.com/svg.image?S.L.Ratio_%7Bsample%7D) can be used as a general quantification of ctDNA for comparing plasma cfDNA samples. This ratio increases when a sample contains the higher contribution of ctDNA.

**Genome-wide Short-fragmented Ratio**

<img src="https://github.com/Pitithat-pu/cfdnakit/wiki/images/wiki/plasma-01-01_OE0290-PED_2LB-055.SLRatio.png" alt="genomewide short-fragment ratio 2LB055" width="800"/>

Normalized

## Package Prerequisites and Installation

Package was tested on R environment 4.0.0. Now only possible way to install this package is via this github repository.

Install prerequisites packages

    if(! "devtools" %in% rownames(installed.packages()))
        install.packages("devtools")
    if(! "BiocManager" %in% rownames(installed.packages()))
        install.packages("BiocManager")

Install cfdnakit package

    library(devtools)  ### use devtools
    install_github("Pitithat-pu/cfdnakit") ### install cfDNAKit 

The installation should work fine without non-zero exit status. Try load cfdnakit package into current R session

    library(cfdnakit) ### Load cfdnakit package

## Usage

Please follow the instructions on [GitHub Wiki page](https://github.com/Pitithat-pu/cfdnakit/wiki)

## Contact

If you have any questions or feedback, please contact us at: Email: [p.puranachot\@dkfz-heidelberg.de](mailto:p.puranachot@dkfz-heidelberg.de); [b.brors\@dkfz-heidelberg.de](mailto:b.brors@dkfz-heidelberg.de)
