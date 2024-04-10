# softmasking of reference sequences by RepeatModeler & RepeatMasker
BuildDatabase -name DATABASE Digitaria_radicosa_v1.fasta
RepeatModeler -database DATABASE -threads 12
RepeatMasker -pa 12 -html -gff -xsmall -nolow -lib DATABASE-families.fa Digitaria_radicosa_v1.fasta
mv Digitaria_radicosa_v1.fasta.masked Digitaria_radicosa_v1.softmasked.fasta

## mapping publicly available RNAseq data
# QC of raw reads
fastp\
    -i read1.fastq.gz -I read2.fastq.gz\
    -o read1_filtered.fastq.gz -O read2_filtered.fastq.gz\
    -h report.html -j report.json -3 -5 -q 20 -n 10 -t 1 -T 1 -l 20 -w 8
# mapping by STAR (perform this for all available samples)
mkdir STAR_index
STAR --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles Digitaria_radicosa_v1.fasta \
    --limitGenomeGenerateRAM 15000000000
STAR --runThreadN 12 \
    --genomeDir STAR_index \
    --readFilesIn read1_filtered.fastq.gz read2_filtered.fastq.gz \
    --readFilesCommand gunzip -c \
    --genomeLoad NoSharedMemory \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix sample

# Gene prediction by BRAKER3
braker.pl \
    --threads=32 \
    --genome=Digitaria_radicosa_v1.softmasked.fasta \
    --bam=Public_RNAseq.sort.bam \
    --prot_seq=merged_aa.fasta

# VUSCO analysis for protein sequences
busco -m protein -i braker/braker.aa -o braker_aa_poales -l poales_odb10 -c 24