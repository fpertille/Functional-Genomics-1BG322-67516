################################################################################################
##### This is the command we used for getting the concordant alignments from the SAM files #####
################################################################################################
module load bioinfo-tools
module load samtools

a=/path/to/bam/list/bam.list.txt
b=/path/to/results
c=name_for_file_containing_concordant_alignments.txt
d=name_for_template_size.txt

for i in $(cat $a);
do
        samtools view $i | grep YT:Z:CP >> $b/$c
done


#################################################################################################
##### This is the command we used for getting the observed template size from the SAM files #####
#################################################################################################

cut -f 9 $b/$c > $b/$d
