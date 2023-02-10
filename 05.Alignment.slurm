###################################################
##### This is the command we used for aligning#####
###################################################

module load bioinfo-tools
module load bowtie2/2.4.5
module load samtools/1.14

a=/path/to/list/of/fastq/
b=/path/to/reference/genome/reference.genome
c=/path/to/the/output/

#Changing the current directory: The script first changes the current working directory to the directory specified by the $a variable.
cd $a
#The script then loops over all files in the current directory that have the .1.fq.gz extension and have names ending with .1.fq.gz. The ls command is used to list all files with the specified extension, and the cut command is used to extract the portion of the filename before the first dot. For each file, the bowtie2 command is executed, aligning the paired-end fastq files specified by the -1 and -2 options to the reference genome specified by the -x option. The -U options specify single-end fastq files, and the -S option specifies the output file in SAM format. The output file is stored in the directory specified by the $c variable, with a name consisting of the extracted portion of the input filename followed by .sam. The --threads option is used to specify the number of threads to use, and the --very-sensitive-local option is used to set the alignment sensitivity.
for i in `ls *.1.fq.gz | cut -d "." -f 1` ;
do bowtie2 --threads 10  --very-sensitive-local -x $b -1 $i\.1.fq.gz -2 $i\.2.fq.gz -U $i\.rem.1.fq.gz -U $i\.rem.2.fq.gz -S $c$i\.sam;
done;


########
#SAMBAM#
########The provided shell script contains four commands, each using the samtools tool to manipulate SAM (Sequence Alignment/Map) format files. The script performs the following operations:##

#Changing the current directory: The script first changes the current working directory to the directory specified by the $a variable.
cd $c

#Converting SAM files to BAM files: The script loops over all files in the current directory that have the .sam extension, and for each file it runs the samtools view command. The -S option is used to indicate that the input is in SAM format, while the -b option is used to output the file in BAM format. The output is saved to a new file with the same name as the input file but with the .sam extension replaced by .bam.
for file in *.sam; do
    samtools view -S -b $file > ${file/%sam/bam}
done

#removing sam files
rm *.sam
#Sorting the BAM files: The script then loops over all BAM files in the current directory and runs the samtools sort command on each one. The -m option is used to specify the maximum memory to be used for sorting, and the sorted output is saved to a new file with the same name as the input file but with .bam replaced by .sorted.bam. 
for file in *.bam; do
    samtools sort -m 768M $file > ${file/%bam/sorted.bam}
done
#The original BAM files are then removed using the ls and xargs rm commands.
ls *.bam| grep -v .sorted.bam$| xargs rm

#Indexing the sorted BAM files: The script then loops over all sorted BAM files in the current directory and runs the samtools index command on each one. The $i".bai" argument specifies the name of the index file to be generated.
for i in *sorted.bam
do
echo "Indexing: "$i
samtools index $i $i".bai"
done
#Finally, the script loops over all sorted BAM files in the current directory and runs the samtools depth command on each one. The output is passed through the awk command, which calculates the average coverage and the total coverage for each file. The average coverage is obtained by dividing the sum of the depths by the number of positions, and the total coverage is obtained by summing the depths.
for file in *.sorted.bam; do
  samtools depth $file |awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}'
done

echo fineshed
