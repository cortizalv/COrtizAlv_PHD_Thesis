#!/bin/bash

star_results=/path/to/star/results
genome_data=/path/to/genome/data
module load star/version 

STAR \
	--runMode genomeGenerate \
	--genomeDir $star_results \
	--genomeFastaFiles $genome_data/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
	--sjdbGTFfile $genome_data/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
	--sjdbOverhang 99 \
	--genomeSAindexNbases 12
