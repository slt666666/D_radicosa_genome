#!/bin/bash

export PATH=${HOME}/local/bin/bwa-0.7.17:${PATH}
export PATH=${HOME}/local/bin:${PATH}

cd ../020-bam

REF_BASE=Dra_ref
if [ ! -e ${REF_BASE}.sa ]; then
	bwa index -p ${REF_BASE} -a is ${REF_BASE}.fa
fi

I=1
for SAMPLE in DG_{1..10} DG_{a..k}; do
	I=$((I+1))
	FASTQ1=../010-fastq/${SAMPLE}.qc.1.fq
	FASTQ2=../010-fastq/${SAMPLE}.qc.2.fq
	BAM=${SAMPLE}.raw.bam
	if [ ! -e ${BAM} ]; then
		bwa mem -t 5 -R "@RG\tID:${I}\tSM:${SAMPLE}" \
			${REF_BASE} ${FASTQ1} ${FASTQ2} > ${BAM} &
	fi
done

wait

