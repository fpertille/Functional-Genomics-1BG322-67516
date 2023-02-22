################################################################################################
##### This is the command we used for getting the concordant alignments from the SAM files #####
################################################################################################
module load bioinfo-tools
module load samtools

a=/path/to/bam
ls -d "$a"/*.bam > $a/bam.list.txt
e=$a/bam.list.txt
b=/path/to/results
c=name_for_file_containing_concordant_alignments.txt
d=name_for_template_size.txt

for i in $(cat $e);
do
        samtools view $i | grep YT:Z:CP >> $b/$c
done

#This code is a shell script that uses a for loop to iterate through each line of a file specified by the variable $a.
#For each line, the script runs a command using the samtools view command to view the contents of a SAM/BAM format file specified by the variable $i.
#The output of the samtools view command is then piped (|) into the grep command to search for a specific pattern, in this case, any line that contains the YT:Z:CP tag. 
#The >> operator appends the output of this command to a file specified by the variable $b/$c.
#Overall, this script appears to be filtering a SAM/BAM file to extract only the lines that contain the YT:Z:CP tag, and saving the output to a specified file location.

#################################################################################################
##### This is the command we used for getting the observed template size from the SAM files #####
#################################################################################################

cut -f 9 $b/$c > $b/$d
rm $c
#This code is a shell command that uses the cut command to extract the ninth field from a file specified by the variable $b/$c, which is the same file that was created in the previous code block.
#The extracted field is then redirected using the > operator to a new file specified by the variable $b/$d.
#Overall, this command appears to be extracting the ninth field (which likely corresponds to a sequence read) from a file, and saving the output to a new file.
