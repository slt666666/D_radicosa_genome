# QC for available NGS reads of Digitaria species
fastp\
    -i read1.fastq.gz -I read2.fastq.gz\
    -o read1_filtered.fastq.gz -O read2_filtered.fastq.gz\
    -h report.html -j report.json -3 -5 -q 20 -n 10 -t 1 -T 1 -l 20 -w 8

# mapping DNA_seq by bwa-mem2 (perform this for all DNA-seq data of Digitaria species)
bwa-mem2 index Digitaria_radicosa_v1.fasta
bwa-mem2 mem Digitaria_radicosa_v1.fasta \
             read1_filtered.fastq.gz \
             read2_filtered.fastq.gz \
             -t 8 | \
samtools sort -@ 8 > sample.sort.bam
samtools index sample.sort.bam

# mapping RNA_seq by STAR (perform this for all RNA-seq data of Digitaria species)
mkdir STAR_index
STAR --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles Digitaria_radicosa_v1.fasta \
    --limitGenomeGenerateRAM 15000000000
STAR --runThreadN 8 \
    --genomeDir STAR_index \
    --readFilesIn read1_filtered.fastq.gz read2_filtered.fastq.gz \
    --readFilesCommand gunzip -c \
    --genomeLoad NoSharedMemory \
    --outFileNamePrefix sample
samtools sort -@ 8 sample.sam > sample.sort.bam
samtools index sample.sort.bam

# Variant calling by GATK
# make GATK index
REF=Digitaria_radicosa_v1.fasta
samtools faidx ${REF}
gatk --java-options "-Xmx16G" CreateSequenceDictionary -R ${REF} -O ${REF%.fasta}.dict

## BQSR for DNAseq (perform this for all DNA-seq data of Digitaria species)
# duplication markup
gatk MarkDuplicates -I sample.sort.bam \
                    -O sample_markdup.bam \
                    -M sample_markdup_metrics.txt
samtools index sample_markdup.bam
gatk AddOrReplaceReadGroups I=sample_markdup.bam O=sample_markdup_addRG.bam \
                            RGLB=4 RGPL=lib1 RGPU=ILLUMINA RGSM=sample
samtools index sample_markdup_addRG.bam
#BQSR
gatk HaplotypeCaller -R ${REF} --emit-ref-confidence GVCF \
                     -I sample_markdup_addRG.bam -O sample.g.vcf

## BQSR for RNAseq (perform this for all RNA-seq data of Digitaria species)
# duplication markup
gatk MarkDuplicates -I sample.sort.bam \
                    -O sample_markdup.bam \
                    -M sample_markdup_metrics.txt
samtools index sample_markdup.bam
gatk AddOrReplaceReadGroups I=sample_markdup.bam O=sample_markdup_addRG.bam \
                            RGLB=4 RGPL=lib1 RGPU=ILLUMINA RGSM=sample
samtools index sample_markdup_addRG.bam
# Splicing handling
gatk SplitNCigarReads -R ../${REF} \
                      -I sample_markdup_addRG.bam -O sample.splitreads.bam
#BQSR
gatk HaplotypeCaller -R ../${REF} --emit-ref-confidence GVCF \
                     -I sample.splitreads.bam -O sample.g.vcf

## build an database and merge gVCF files
gvcf_files=""
for gvcf_file in `ls *.g.vcf`; do
 gvcf_files=${gvcf_files}"-V ${gvcf_file} "
done
echo -e "chr01\nchr02\nchr03\nchr04\nchr05\nchr06\nchr07\nchr08\nchr09" > intervals.list
gatk GenomicsDBImport -R ${REF} ${gvcf_files} -L intervals.list \
                      --genomicsdb-workspace-path gvcfs_db
gatk GenotypeGVCFs -R ${REF} -V gendb://gvcfs_db -O merged.vcf
gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include SNP -O merged_snps.vcf
gatk VariantFiltration -R ${REF} -V merged_snps.vcf -O merged_snps_filtered.vcf \
                      -filter "QD < 2.0" --filter-name "QD2"       \
                      -filter "QUAL < 30.0" --filter-name "QUAL30" \
                      -filter "SOR > 4.0" --filter-name "SOR4"     \
                      -filter "FS > 60.0" --filter-name "FS60"     \
                      -filter "MQ < 40.0" --filter-name "MQ40"     \
                      -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                      -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include INDEL -O merged_indels.vcf
gatk VariantFiltration -R ${REF} -V merged_indels.vcf -O merged_indels_filtered.vcf \
                      -filter "QD < 2.0" --filter-name "QD2"       \
                      -filter "QUAL < 30.0" --filter-name "QUAL30" \
                      -filter "FS > 200.0" --filter-name "FS200"   \
                      -filter "SOR > 10.0" -filter-name "SOR10"    \
                      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

# Apply BQSR for DNAseq (perform this for all DNA-seq data of Digitaria species)
gatk BaseRecalibrator -R ${REF} -I sample_markdup_addRG.bam \
                    --known-sites ../bqsr/merged_snps_filtered.vcf \
                    --known-sites ../bqsr/merged_indels_filtered.vcf \
                    -O sample_recal_data.table
gatk ApplyBQSR -R ${REF} -I sample_markdup_addRG.bam -bqsr sample_recal_data.table -O sample_bqsr.bam

# Variant call for DNAseq (perform this for all DNA-seq data of Digitaria species)
gatk HaplotypeCaller -R ${REF} --emit-ref-confidence GVCF --native-pair-hmm-threads 12 \
                     -I sample_bqsr.bam -O sample.g.vcf

# Apply BQSR for RNAseq (perform this for all RNA-seq data of Digitaria species)
gatk BaseRecalibrator -R ${REF} -I sample.splitreads.bam \
                    --known-sites ../bqsr/merged_snps_filtered.vcf \
                    --known-sites ../bqsr/merged_indels_filtered.vcf \
                    -O sample_recal_data.table
gatk ApplyBQSR -R ${REF} -I sample.splitreads.bam -bqsr sample_recal_data.table -O sample_bqsr.bam

# Variant call for RNAseq (perform this for all RNA-seq data of Digitaria species)
gatk HaplotypeCaller -R ${REF} --emit-ref-confidence GVCF --native-pair-hmm-threads 12 \
                     -I sample_bqsr.bam -O sample.g.vcf

# list up all gVCF files
gvcf_files=""
for gvcf_file in `ls *.g.vcf`; do
 gvcf_files=${gvcf_files}"-V ${gvcf_file} "
done
# Variant call using all samples
echo -e "chr01\nchr02\nchr03\nchr04\nchr05\nchr06\nchr07\nchr08\nchr09" > intervals.list
gatk GenomicsDBImport -R ${REF} ${gvcf_files} -L intervals.list \
                     --genomicsdb-workspace-path gvcfs_db
gatk GenotypeGVCFs -R ${REF} -V gendb://gvcfs_db -O merged.vcf
gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include SNP -O merged_snps.vcf
gatk VariantFiltration -R ${REF} -V merged_snps.vcf -O merged_snps_filtered.vcf \
                      -filter "QD < 2.0" --filter-name "QD2"       \
                      -filter "QUAL < 30.0" --filter-name "QUAL30" \
                      -filter "SOR > 4.0" --filter-name "SOR4"     \
                      -filter "FS > 60.0" --filter-name "FS60"     \
                      -filter "MQ < 40.0" --filter-name "MQ40"     \
                      -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                      -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include INDEL -O merged_indels.vcf
gatk VariantFiltration -R ${REF} -V merged_indels.vcf -O merged_indels_filtered.vcf \
                      -filter "QD < 2.0" --filter-name "QD2"       \
                      -filter "QUAL < 30.0" --filter-name "QUAL30" \
                      -filter "FS > 200.0" --filter-name "FS200"   \
                      -filter "SOR > 10.0" -filter-name "SOR10"    \
                      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

# Select biallelic variants
gatk SelectVariants \
   -R ${REF} \
   -V merged_snps_filtered.vcf  \
   -O merged_snps_filtered_biallelic.vcf  \
   --restrict-alleles-to BIALLELIC \
   --exclude-filtered  \

# extract variants with deapth > 1
bcftools view -i 'MIN(FMT/DP)>1' merged_snps_filtered_biallelic.vcf > Genetic_variants_Digitaria_species.vcf
