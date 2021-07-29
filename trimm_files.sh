#!/bin/bash

data_path=/path/to/data #path where the .fq data is stored 
names_path=/path/to/a/file/with/sample/names #path to a file that contains the names of the samples 
trimm_results=/path/to/trimmomatic/results #path where the trimmomatic output will be stored

find $data_path/*_1.fq |
while read fname 
 do
   base=$(basename ${fname} _1.fq)
   java -Xmx80g -jar /cm/shared/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 10 -phred33 $fname $data_path/${base}_2.fq \
                $trimm_results/${base}_1.trim.fastq $trimm_results/${base}_1.trim.unpaired.fastq \
                $trimm_results/${base}_2.trim.fastq $trimm_results/${base}_2.trim.unpaired.fastq \
                ILLUMINACLIP:$names_path/TruSeq3-PE-2.fa:2:15:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:30
 done
