# Scripts for Puerto Rico Honey bee brain RNASeq Analysis

This repository includes scripts used to process and analyze RNASeq data from honey bee brains.

(C) Carlos A. Ortiz Alvarado 2021

## Scripts 

1. `fastqc_files.sh`: This script processes raw FASTQ reads using *fastqc*
2. `multiqc_files.sh`: This script compiles outputs from *fastqc* into a single quality report using *multiqc*
3. `trimm_files.sh`: This script removes low-quality reads from the raw FASTQ data using *trimmomatic*
4. `index_star.sh`: This script generates a genome and annotation index using *star*
5. `align_star_untrimmed.sh`: This script aligns FASTQ reads of all samples  and generates a count data matrix using *star*
