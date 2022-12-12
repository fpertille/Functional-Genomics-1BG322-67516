#! /bin/bash

##################################################################################
########### In this first part the GFF file is going to be created ###############
##################################################################################
# First the GFF file has to be converted to unix format
dos2unix <input>

# The gff file has a header that need to be extracted before doing the first two chunks of commands

# We need to put a unique identifier for each of the windows
# Creating a file with the ID name for as many times as cut you have in the ggf file
START=1
END=$(wc -l < ./in_silico_windows_mm39.gff)

for (( c=$START; c<=$END; c++ ));
  do
    echo "ID=$c;gene_id=$c;"
  done > columns1.txt

# Now we just need to paste it along with the rest of the columns
paste in_silico_windows_<input>.gff columns1.txt > in_silico_windows.curated.<input>.gff
# Add the header at the beginning of the GFF file

##################################################################################
############### In this second part the count matrix is created ##################
##################################################################################

# FeatureCounts recognices the windows in the BAM file and produces a count matrix, this count matrix is going to be later used to do the statistic to see the differential methylated windows
featureCounts -t sequence_feature -g gene_id -p -B -a /<path_GFF_file>/<GFF_file> -o /<path_output>/counts.<output>.txt /<path_input>/*.bam

