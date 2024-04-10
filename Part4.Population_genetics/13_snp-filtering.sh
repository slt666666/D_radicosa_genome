#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}
cd ../030-vcf/

DCI_BCF=Dci2.bcf
DCI_FLT_VCF=Dci2.flt.vcf
bcftools view --with-header ${DCI_BCF} |
vcftools --vcf - --maf 0.05 --minGQ 20 --minDP 5 --maxDP 50 \
	--max-alleles 2 --recode --stdout |
vcftools --vcf - --max-missing 0.8 --recode --stdout |
bcftools +prune -n 1 -w 1kb --random-seed 543 -Ov - |
bcftools annotate \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	> ${DCI_FLT_VCF} &

ALL_BCF=DG2.bcf
ALL_FLT_VCF=DG2.flt.vcf
bcftools view --with-header ${ALL_BCF} |
vcftools --vcf - --maf 0.05 --minGQ 20 --minDP 5 --maxDP 50 \
	--max-alleles 2 --recode --stdout |
vcftools --vcf - --max-missing 0.8 --recode --stdout |
bcftools +prune -n 1 -w 1kb --random-seed 543 -Ov - |
bcftools annotate \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	> ${ALL_FLT_VCF} &

wait

mkdir -p Dci2 DG2
mv ${DCI_FLT_VCF} Dci2/
mv ${ALL_FLT_VCF} DG2/
