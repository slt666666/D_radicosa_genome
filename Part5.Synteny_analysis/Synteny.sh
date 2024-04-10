# Synteny analysis by MCScanX
makeblastdb -in Rice_reference/IRGSP-1.0_protein_2024-01-11.fasta -out Rice_DB -dbtype prot -parse_seqids
makeblastdb -in Fonio_reference/protein.fasta -out Fonio_DB -dbtype prot -parse_seqids
makeblastdb -in Digitaria_radicosa_v1.protein.fasta -out Radicosa_DB -dbtype prot -parse_seqids

blastp -query Rice_reference/IRGSP-1.0_protein_2024-01-11.fasta -db Radicosa_DB -evalue 1e-10 -num_alignments 5 -num_descriptions 5 -outfmt 6 -num_threads 18 -out Radicosa_vs_Rice.blast2
blastp -query Fonio_reference/protein.fasta -db Radicosa_DB -evalue 1e-10 -num_alignments 5 -num_descriptions 5 -outfmt 6 -num_threads 18 -out Radicosa_vs_Fonio.blast2
blastp -query Digitaria_radicosa_v1.protein.fasta -db Fonio_DB -evalue 1e-10 -num_alignments 5 -num_descriptions 5 -outfmt 6 -num_threads 18 -out Radicosa_vs_Fonio.blast1
blastp -query Digitaria_radicosa_v1.protein.fasta -db Rice_DB -evalue 1e-10 -num_alignments 5 -num_descriptions 5 -outfmt 6 -num_threads 18 -out Radicosa_vs_Rice.blast1

cat radicosa.gff fonio.gff > Radicosa_vs_Fonio.gff
cat radicosa.gff rice.gff > Radicosa_vs_Rice.gff

cat Radicosa_vs_Rice.blast1 Radicosa_vs_Rice.blast2 > Radicosa_vs_Rice.blast
cat Radicosa_vs_Fonio.blast1 Radicosa_vs_Fonio.blast2 > Radicosa_vs_Fonio.blast

MCScanX ./Radicosa_vs_Fonio -s 10 -b 2
MCScanX ./Radicosa_vs_Rice -s 10 -b 2