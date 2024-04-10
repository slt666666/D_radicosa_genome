#!/bin/bash

# Call SNPs again without samples presumably atttibuted to
# wrong species (DG_a & G_c)

export PATH=${HOME}/local/bin:${PATH}

cd ../030-vcf/
REFERENCE=Dra_ref.fa
BCF_DCI=Dci2.bcf
BCF_ALL=DG2.bcf

bcftools mpileup \
	-a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/SCR \
	-f ${REFERENCE} -d 10000 -Ov \
	--samples `echo DG_b DG_{d..k} | tr ' ' ','` \
	DG.bam |
bcftools call -m --annotate GQ,GP -Ov - |
bcftools view -v snps -Ob - \
> ${BCF_DCI} &

bcftools mpileup \
	-a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/SCR \
	-f ${REFERENCE} -d 10000 -Ov \
	--samples `echo DG_{1..10} DG_b DG_{d..k} | tr ' ' ','` \
	DG.bam |
bcftools call \
	-m --annotate GQ,GP -Ov - |
bcftools view -v snps -Ob - \
> ${BCF_ALL} &
wait

