#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}

cd ../010-fastq
mkdir -p fastqc

for ACC in DG_{1..10} DG_{a..k}; do
	fastqc -o fastqc -f fastq -t 6 ${ACC}.qc.{1,2}.fq &
done
wait

