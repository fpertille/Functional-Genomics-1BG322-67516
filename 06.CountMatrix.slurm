#! /bin/bash

##################################################################################
######################### The count matrix is created ############################
##################################################################################
#path to gff file
i=path_to_gff_file
# path to output folder
j=path_to_output
#path to the input folder (BAM files)
k=path_to_input
# name of the reference genome to make it unique
b=reference_genome_name (like mm10)
# name of the gff file
h=file.curated.gff

module load bioinfo-tools
module load subread

# FeatureCounts recognices the windows in the BAM file and produces a count matrix, this count matrix is going to be later used to do the statistic to see the differential methylated windows
featureCounts -t sequence_feature -g gene_id -p -B -a $i/$h -o $j/counts.$b.txt $k/*.bam

# featureCounts tool to quantify the reads mapped to features, such as genes or exons, in BAM files specified by $k/*.bam. 
#The -t option specifies the type of features, 
#the -g option specifies the identifier for the features, 
#the -p option enables the counting of paired-end reads, 
#the -B option enables the counting of both reads for paired-end reads that have one end mapped, 
#and the -a option specifies the annotation file in the GTF format. 
#The -o option specifies the output file. 
#The counts for each feature are stored in the specified output file in the tab-delimited format.

