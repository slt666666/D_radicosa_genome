#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}

cd ../010-fastq
mkdir -p fastqc

for FILE in DG_{1..10}.{1,2}.fq.gz; do
	gzip -df ${FILE} &
done
wait

for FILE in DG_{a..k}.{1,2}.fq.gz; do
	gzip -df ${FILE} &
done
wait

for ACC in DG_{1..10}; do
	fastqc -o fastqc -f fastq -t 6 ${ACC}.{1,2}.fq &
done
wait

for ACC in DG_{a..k}; do
	fastqc -o fastqc -f fastq -t 6 ${ACC}.{1,2}.fq &
done
wait

