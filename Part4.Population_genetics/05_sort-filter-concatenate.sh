#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}

cd ../020-bam/

for SAMPLE in DG_{1..10} DG_{a..k}; do
	RAW_BAM=${SAMPLE}.raw.bam
	FLT_BAM=${SAMPLE}.flt.bam
	if [ ! -e ${FLT_BAM} ]; then
		samtools sort --threads 1 ${RAW_BAM} |
		samtools view --threads 1 --with-header \
			 --expr 'flag.paired == 1 && flag.proper_pair == 2' - |
		samtools view --threads 1 --with-header --min-MQ 40 - |
		samtools view --threads 1 --with-header - |
			grep -v "SA:Z:" | grep -v "XA:Z:" |
		samtools view --threads 1 --bam --with-header - > ${FLT_BAM} &
	fi
done
wait

ALL_BAM=DG.bam
samtools merge ${ALL_BAM} DG_{1..10}.flt.bam DG_{a..k}.flt.bam

