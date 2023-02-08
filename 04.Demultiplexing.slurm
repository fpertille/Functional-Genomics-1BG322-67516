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

process_radtags -P -1 $h/$g -2 $h/$f -o $a -b $c -e pstI -r -q -c --len-limit 20
