# cfdnakit : an R package for fragmentation analysis of cfDNA and copy-number alteration calling

## Introduction

This package provides mostly basic functions to analysis cell-free DNA (cfDNA) with next-generation sequencing. The package focuses on extracting length of cfDNA,  and genome-wide copy-number alteration estimated by the short-fragmented cfDNA using shallow whole-genome sequencing data (~0.3X or more).

## Features

The scope of this R package is to analyse the length of cfDNA fragments. The package simplifies the process of extracting length of fragments from a BAM file and provides basic functions to explore this characteristic of cfDNA with low-coverage whole-genome sequencing data. Moreover, this package utilizies the quantity of short-fragmented cfDNA to infer copy-number alterations and estimate the percentage of tumor-derived cfDNA.

### Fragment length distribution and comparison

Package provides a single function to extract fragment length of cfDNA in the sample. Making a fragment-length distribution plot of multiple samples is easy. cfdnakit also extracted the short-fragment ratio representing the amount of short-fragmented cfdNA in the sample. It can be used for comparison between groups of sample (e.g. healthy vs patient, stable vs advance stage) or for quality control inspection. 

This plot shows the fragment-length ditribution of cfDNA from a PDX sample. One line represent tumor-derived (human) cfDNA , another represents non-tumor (mouse) cfDNA.

<img src="https://github.com/Pitithat-pu/cfdnakit/wiki/images/wiki/fragment_xenograft_cfdnakit.png" title="cfdnakit fragment length distribution" alt="fragment_length_distribution_cfdnakit" width="550"/>

### Copy-number calling and tumor fraction estimation from short-fragmented cfDNA

The figure below shows the overview of the CNV-calling procedure. The amount of short-fragmented cfDNA per windows are normalised and compared to a Panel-of-Normal (control). Segmentation is performed using the PSCBS package. A CPA score (Raman, Lennart, et al. 2020) were calculated to estimate the copy number tumor burden.

<img src="https://raw.githubusercontent.com/wiki/Pitithat-pu/cfdnakit/images/cfdnakit_workflow.png" title="cfdnakit cnv calling workflow" alt="cfdnakit_workflow" width="600"/>

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
