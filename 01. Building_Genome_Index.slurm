#Violeta, can you fix it the way you did with the otherones?
bowtie2-build /proj/gbs_medip/forgeGenome/Galgal7/ncbi_dataset/data/GCF_016699485.2/seqs_srcdir/newgenomeclean.fa /proj/gbs_medip/uppstore2017266/analysis/genomes/Galgal7/Galgal7

#The command bowtie2-build is used to create an index of a reference genome in order to align reads from sequencing data to the genome. 
#The command takes two arguments: the path to the reference genome file in FASTA format, and the name of the output index.
#In this case, the reference genome file is located at /proj/gbs_medip/forgeGenome/Galgal7/ncbi_dataset/data/GCF_016699485.2/seqs_srcdir/newgenomeclean.fa and the output index will be saved as /proj/gbs_medip/uppstore2017266/analysis/genomes/Galgal7/Galgal7.
#So, the command is generating a Bowtie2 index for the genome assembly of the chicken (Gallus gallus) with the version GCF_016699485.2. The output index file will be saved under /proj/gbs_medip/uppstore2017266/analysis/genomes/Galgal7/ with the prefix Galgal7. This index will be used later on to align sequencing reads from this organism to the reference genome.
