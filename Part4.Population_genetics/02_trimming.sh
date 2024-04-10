#!/bin/bash

export PATH=${HOME}/local/bin:${PATH}

# adapter.fasta was created from the sequences of
# the MIG-seq 2nd primers and the indeces.

cd ../010-fastq/
mkdir -p fastp

for ID in {1..10} {a..k}; do
fastp --in1 DG_${ID}.1.fq --out1 DG_${ID}.qc.1.fq \
	--in2 DG_${ID}.2.fq --out2 DG_${ID}.qc.2.fq \
	--adapter_fasta adapters.fasta \
	--trim_front1 17 --trim_tail1 5 \
	--trim_front2 17 --trim_tail2 5 \
	--trim_poly_x \
	--cut_front --cut_tail --cut_mean_quality 30 \
	--average_qual 30 \
	--length_required 20 \
	--json fastp/DG_${ID}.fastp.json \
	--html fastp/DG_${ID}.fastp.html \
	--thread 4 &
done

