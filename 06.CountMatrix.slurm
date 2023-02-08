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

