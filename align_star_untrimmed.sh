#!/bin/bash

star_results=/path/to/star/results 
data_path=/path/to/data
genome_data=/path/to/genome/data
align_output=/path/to/alignment/results
module load star/version_number

find $data_path/*_1.fq |
while read fname 
 do
	base=$(basename ${fname} _1.fq)
	STAR \
	--runMode alignReads \
	--runThreadN 12 \
	--genomeDir $star_results \
	--readFilesIn $data_path/${base}_1.fq $data_path/${base}_2.fq \
	--sjdbGTFfile $genome_data/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
	--sjdbOverhang 99 \
	--outFileNamePrefix $align_output/${base} \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outTmpDir ${SLURM_JOB_ID}\
done
