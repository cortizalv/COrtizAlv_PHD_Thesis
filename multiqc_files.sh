#!/bin/bash

fastqc_dir=/path/to/fastqc/results #path where the fastqc data is stored 
out_dir=/path/to/output/directory #path where the multiqc output will be stored 

module load multiqc/version_number 

multiqc $fastqc_dir/*_fastqc.zip -o $out_dir
