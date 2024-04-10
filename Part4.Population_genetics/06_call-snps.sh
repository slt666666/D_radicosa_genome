#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}

cd ../030-vcf/
REFERENCE=Dra_ref.fa
BCF_DRA=Dra.bcf
BCF_DCI=Dci.bcf
BCF_ALL=DG.bcf

bcftools mpileup \
	-a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/SCR \
	-f ${REFERENCE} -d 10000 -Ov \
	--samples `echo DG_{1..10} | tr ' ' ','` \
	DG.bam |
bcftools call -m --annotate GQ,GP -Ov - |
bcftools view -v snps -Ob - \
> ${BCF_DRA} &

bcftools mpileup \
	-a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/SCR \
	-f ${REFERENCE} -d 10000 -Ov \
	--samples `echo DG_{a..k} | tr ' ' ','` \
	DG.bam |
bcftools call -m --annotate GQ,GP -Ov - |
bcftools view -v snps -Ob - \
> ${BCF_DCI} &

bcftools mpileup \
	-a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/SCR \
	-f ${REFERENCE} -d 10000 -Ov \
	DG.bam |
bcftools call \
	-m --annotate GQ,GP -Ov - |
bcftools view -v snps -Ob - \
> ${BCF_ALL} &
wait

