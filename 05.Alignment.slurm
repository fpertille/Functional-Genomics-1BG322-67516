###################################################
##### This is the command we used for aligning#####
###################################################

module load bioinfo-tools
module load bowtie2/2.4.5

a=/path/to/list/of/fastq/fastq.list.txt
b=/path/to/reference/genome/reference.genome.fa

for i in $(cat $a) ;

do bowtie2 --threads 10  --very-sensitive-local -x $b -1 $i -2 ${i%%.1.fq.gz}.2.fq.gz -U ${i%%.1.fq.gz}.rem.1.fq.gz -U ${i%%.1.fq.gz}.rem.2.fq.gz -S ${i%%.1.fq.gz}.sam ;

done
