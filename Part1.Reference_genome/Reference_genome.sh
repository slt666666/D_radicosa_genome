# QC for illumina short reads
fastp\
 -i Digitaria_leaf_R1.fastq.gz -I Digitaria_leaf_R2.fastq.gz\
 -o Digitaria_leaf_R1_filtered.fastq.gz -O Digitaria_leaf_R2_filtered.fastq.gz\
 -h report.html -j report.json -3 -5 -q 30 -u 60 -e 20 -n 10 -t 1 -T 1 -l 20 -w 8

# Hifiasm
hifiasm -o Digitaria.asm -t 24 Digitaria.fastq.gz 2> hifiasm.log
awk '/^S/{print ">"$2;print $3}' Digitaria.asm.bp.p_ctg.gfa > Digitaria.asm.p_ctg.fa

# mapping HiFi reads to hifiasm output
meryl count k=15 output merylDB Digitaria.asm.p_ctg.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -t 24 -W repetitive_k15.txt -ax map-pb Digitaria.asm.p_ctg.fa Digitaria.fastq.gz|samtools sort -o hifi_sort.bam
samtools index -@ 24 hifi_sort.bam

## Polishing by nextPolish2
# calculate k-mer
# 21-mer
yak count -o k21.yak -k 21 -b 37 <(zcat Illumina_short_read/*filtered.fastq.gz) <(zcat Illumina_short_read/*filtered.fastq.gz)
# 31-mer
yak count -o k31.yak -k 31 -b 37 <(zcat Illumina_short_read/*filtered.fastq.gz) <(zcat Illumina_short_read/*filtered.fastq.gz)
# nextPolish2
nextPolish2 -t 8 hifi_sort.bam .Digitaria.asm.p_ctg.fa k21.yak k31.yak > Digitaria_radicosa.hifiasm_p_ctg.nextPolish2.fasta

# Purge duplicates
python run_purge_dups.py -p bash config.json ~/purge_dups/src Digitaria_radicosa

# Reference guided assembly by RagTag -> change output name to Digitaria_radicosa_v1.fasta
ragtag.py scaffold -o ragtag_fonio_B_genome Digitaria_exilis.Fonio_CM05836.dna.toplevel_B_genome.fa Digitaria_radicosa.hifiasm_p_ctg.nextPolish2.purged.fa

# BUSCO analysis
busco -m genome -i Digitaria_radicosa_v1.fasta -o ref_genome_poales -l poales_odb10 -c 12 --update-data
