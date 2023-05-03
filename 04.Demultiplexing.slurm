#!/bin/bash -l
#SBATCH -A uppmax2022-2-30
#SBATCH -M snowy
#SBATCH -p core -n 1
#SBATCH --mem=64gb
#SBATCH -t 24:00:00
#SBATCH -J trimm
#SBATCH --error /error_path/01.DemultiplexingLane1.err
#SBATCH --output /output_path/01.DemultiplexingLane1.out

#######################################################################################################################################
##### This is the command we used for demultiplexing, you would need to change the input file as we don't want to make it chaotic #####
#######################################################################################################################################

module load bioinfo-tools
module load Stacks/2.62
#rescued from the 03 script, you should have the same name
h=/path/to/filter/folder
g=MeDIP.R1.nopoliG.fastq.gz #name of the R1 input
f=MeDIP.R2.nopoliG.fastq.gz #name of the R2 input
a=/path/to/demultiplexed/folder
c=/path/to/barcode/file/barcode.txt

process_radtags --threads 4  -P -1 $h/$g -2 $h/$f -o $a -b $c -e pstI -r -q -c --len-limit 20

#--threads 4 indicates that the program will use 4 threads to process the data.
#-P specifies that the input files are paired-end reads.
#-1 $h/$g and -2 $h/$f indicate the paths to the input files for the forward and reverse reads, respectively. The $h variable is a placeholder for the path to the filter folder, and $g and $f are the names of the input files for the forward and reverse reads.
#-o $a specifies the path to the output directory for the processed files. The $a variable is a placeholder for the path to the demultiplexed folder.
#-b $c indicates the path to the file containing the barcode sequences used to demultiplex the reads. The $c variable is a placeholder for the path to the barcode file.
#-e pstI specifies the restriction enzyme used to generate the RAD-seq library.
#-r ,--rescue — rescue barcodes and RAD-Tag cut sites.
#-q specifies that the program should perform quality filtering on the reads.
##by default, process_radtags uses a sliding window quality filter to remove reads with average quality scores below a certain threshold. This threshold is set to 20 by default, but it can be adjusted using the -y flag.
###The sliding window quality filter works by calculating the average quality score of a sliding window of a fixed length (by default, 15 bp) along each read. If the average quality score falls below the specified threshold, the read is discarded.
#-c --clean — clean data, remove any read with an uncalled base.
#--len-limit 20 specifies the maximum length of the reads to be retained after trimming. Any reads longer than 20 bp will be discarded.
