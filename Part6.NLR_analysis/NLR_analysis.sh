# perform InterProscan
interproscan.sh -i Digitaria_radicosa_v1.protein.fasta -f gff3 -o Digitaria_radicosa_v1.protein.interproscan.gff -cpu 16 -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles

# perform NLRtracker to annotate NLRs
i=Digitaria_radicosa_v1.protein
NLRtracker.sh -s ../Step10.rename_genes/${i}.fasta -o ${i}.NLRtracker -i ${i}.interproscan.gff

# alignment NLRs
cat Poales_functionally_validated_NLR.fasta D_radicosa_NLR.fasta > Poales_dedup_Radicosa_NLR.fasta
mafft Poales_dedup_Radicosa_NLR.fasta > Poales_dedup_Radicosa_NLR_alignment.fasta

# manually extracted NBARC domain from alignment results
# -> make phylogenetic tree by RAxML
raxmlHPC-PTHREADS-AVX2 \
-T 12 \
-# 100 \
-f a \
-x 1024 \
-m PROTGAMMAAUTO \
-s Poales_dedup_Radicosa_NLR_alignment_NBARC.fasta \
-p 121 \
-n Poales_dedup_Radicosa_NLR