#! /bin/bash

##################################################################################
######################### The count matrix is created ############################
##################################################################################

module load subread
# FeatureCounts recognices the windows in the BAM file and produces a count matrix, this count matrix is going to be later used to do the statistic to see the differential methylated windows
featureCounts -t sequence_feature -g gene_id -p -B -a /<path_GFF_file>/<GFF_file> -o /<path_output>/counts.<output>.txt /<path_input>/*.bam

