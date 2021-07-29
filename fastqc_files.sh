#!/bin/bash

data_dir=/path/to/data #Directory where the data is found
results_dir=/path/to/results #path where the results of the fastqc will be stored 
module load fastqc/version_number #load the fastqc module

find $data_dir -type f -name "*.fq" | 
while read fname
do
	fastqc ${fname} -o $results_dir
done
