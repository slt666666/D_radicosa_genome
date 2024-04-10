# removing and counting invariant sites from a PHYLIP file to use for the RAxML ascertainment bias correction
python3 vcf2phylip.py --input Genetic_variants_Digitaria_species.vcf
python3 ascbias.py -p Genetic_variants_Digitaria_species.min4.phy

# selecting the best-fit model of evolution for DNA alignments.
modeltest-ng -i out.phy

# make phylogenetic tree by RAxML
/home/b/b34531/t.sakai/tools/standard-RAxML/raxmlHPC-PTHREADS-AVX2 \
-T 12 \
-# 1000 \
-f a \
-x 1024 \
-m GTRGAMMA \
-s out.phy \
-p 121 \
-n out.phy