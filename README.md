# cfdnakit : an R package for fragmentation analysis of cfDNA and copy-number alteration calling

## Introduction

This package provides mostly basic functions to analysis cell-free DNA (cfDNA) with next-generation sequencing. This first version will focus on analyzing shallow whole-genome sequence including fragment length analysis, genome-wide fragmentation, and copy-number alteration estimated by short cell-free DNA fragments.

![cfdnakit_workflow](https://raw.githubusercontent.com/wiki/Pitithat-pu/cfdnakit/images/cfdnakit_workflow.png)

## Features

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
