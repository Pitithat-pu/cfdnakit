# cfDNAKit : an R package for next-generation sequencing analysis focusing on fragmentation aspect of liquid biopsy cfDNA 

## Introduction

This package provides mostly basic functions to analysis cell-free DNA (cfDNA) with next-generation sequencing. This first version will focus on analyzing shallow whole-genome sequence including fragment length analysis, genome-wide fragmentation, and copy-number alteration estimated by short cell-free DNA fragments.

## Package Prerequisites and installation
Package were run and tested on R environment 4.0.0 .
Now only possible way to install this package is via this github repository which require package devtools and BiocManager.

Install prerequisites packages
```
if(! "devtools" %in% rownames(installed.packages()))
	install.packages("devtools")
if(! "BiocManager" %in% rownames(installed.packages()))
	install.packages("BiocManager")
```
Install cfDNAKit package
```
library(devtools)  ### use devtools
install_github("Pitithat-pu/cfdnakit") ### install cfDNAKit 
```
The installation should work fine without non-zero exit status.
## Contact
If you have any questions or feedback, please contact us at:\
Email: p.puranachot@dkfz-heidelberg.de; b.brors@dkfz-heidelberg.de
