##############################################################################################################################################################################
##### This is the command we used for getting the concordant alignments from the SAM files , you would need to change the input file as we don't want to make it chaotic #####
##############################################################################################################################################################################

grep YT:Z:CP C13F2_21_CO_sorted_reheader.sam > concordant_alignment_C13F2_21_CO_sorted_reheader

###############################################################################################################################################################################
##### This is the command we used for getting the observed template size from the SAM files , you would need to change the input file as we don't want to make it chaotic #####
###############################################################################################################################################################################

cut -f 9 concordant_alignment_C13F2_21_CO_sorted_reheader > outputfile
