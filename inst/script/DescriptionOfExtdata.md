---
title: Description of external data in cfdnakit package
date: 28.05.2023
author: Pitithat Puranachot
---

### AnnotationDataFrame_from_QDNAseq_hg19_1000k.rds and AnnotationDataFrame_from_QDNAseq_hg38_1000k.rds
These RDS files are pre-download AnnotatedDataFrame objects obtained from QDNAseq::getBinAnnotations function for hg19 and hg38 reference genome from QDNAseq.
```
QDNAseq::getBinAnnotations(binSize=1000,genome="hg19")
```
They are essential for the analysis when we split fragment-length information into 1 megabase non-overlapping bins.
Other bin sizes, like 500K and 100K, will be downloaded on the fly. However, you can obtain the AnnotatedDataFrame object, save as RDS file named in the same naming format, and put it in the package's extdata directory to reuse the object.

### BH01_CHR15.SampleBam.rds
This file is a SampleBam object of cfDNA of healthy individuals cfdna BH01 published in [Snyder et.al](10.1016/j.cell.2015.11.050).
The total sequence data were aligned with BWA-MEM onto the human reference genome (hg19) and only reads on chromosome 15 are extracted. Then the sequence alignment were downsampled around 0.005 fold to limit the size of the file to comply the Bioconductor file-size limit.
This file were used to demonstrate how to integrate multiple sample into fragment-length distribution plot shown in package's vignette section **Analyse the Fragment Length Distribution**. 

### example_patientcfDNA_SampleBam.RDS
This file stores a SampleBam object created by function read_bamfile from a cfDNA sample of a cancer patient from pediatric cancer cohort at KiTZ - Hopp-Kindertumorzentrum Heidelberg. The alignment file were downsampled to limit the size of the file to comply the Bioconductor file-size limit.
This file were also used to demonstrate how to integrate multiple sample into fragment-length distribution plot shown in package's vignette section **Analyse the Fragment Length Distribution**.

### example_patientcfDNA_SampleFragment.RDS
This file stores a SampleFragment object created by function get_fragment_profile from the same patient of file example_patientcfDNA_SampleBam.RDS but from a higher sequencing coverage data. in order to properly demonstrate the analyse of genome-wide fragment-length profile.
This file is for demonstrating how to extract fragment-length information per non-overlapping bins shown in package's vignette section **Plot Genome-wide Short-fragmented Ratio**. It also used on example documentation of many functions such as plot_sl_ratio, plot_transformed_sl, get_zscore_profile and etc.

### ex.healthy1.rds and ex.healthy2.rds

Both files were used for example documentation of the function create_PoN.R. They are SampleFragment objects from cfDNA sample of two healthy individuals from pediatric cancer cohort at KiTZ - Hopp-Kindertumorzentrum Heidelberg. Their alignment files were read by read_bamfile function, producing SampleBam object, and then by get_fragment_profile function.

### ex.plasma.bam
This sequence alignment file stored sequencing data of a cfDNA sample of a patient from pediatric cancer cohort at KiTZ - Hopp-Kindertumorzentrum Heidelberg. The file were downsampled around 0.1 fold. It is used for function demonstration in package's vignette and example documentation of the function read_bamfile. 

### ex.PoN.rds
This files is an object to store the result dataframe of create_PoN function. This Panel-of-Normal information was created from 10 cfDNAs of healthy individuals from pediatric cancer cohort at KiTZ - Hopp-Kindertumorzentrum Heidelberg. It will be used for providing example document of many packages' functions such as get_zscore_profile, plot_transformed_sl.R, etc. The PoN information is necessary for downstream analysis.

### hg19_centromere.tsv.gz
This file stores regions of centromere of human chromosomes. The information were obtained from [UCSC genome browser](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz) where centromere regions are annotated with "acen" in the last column. So only regions with "acen" are included into this file. This information were used to filter read mapping onto the centromere regions.

### hg19_chrTotalLength.tsv and hg38_chrTotalLength.tsv
They are storing number of base-pair per chromosome of the reference genome hg19 and hg38. They are used as an information for plotting result of genome-wide data analysis such as plotting SL.Ratio per non-overlapping bin. The source is from the NCBI Genome Reference Consortium [hg19](https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37) and [hg38](https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38).

### wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz
This file contains DacMapability regions for the human reference hg19. All reads mapped onto those regions that will be excluded from analysis.
[UCSC genome browser](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeDacMapabilityConsensusExcludable.txt.gz)
