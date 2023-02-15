##################################################################################################################################
##### This is the command we used for filtering, here you can see the filter for trimming adapters and trimming of poliG > 2 #####
##################################################################################################################################

module load bioinfo-tools
module load fastp
i=/path/to/R1_001.fastq.gz
j=/path/to/R2_001.fastq.gz
h=/path/to/filter/folder
g=MeDIP.R1.nopoliG.fastq.gz #name of the R1 output
f=MeDIP.R2.nopoliG.fastq.gz #name of the R2 output
w=MeDIP.nopoliG.html #name of the report

fastp --thread 4 -z 4 --in1 $i --in2 $j --out1 $h/$g --out2 $h/$f -z 4 --detect_adapter_for_pe --poly_g_min_len 2 -h $h/$w


# The -z option is used to specify the number of worker threads, the --detect_adapter_for_pe option is used to enable adapter detection for paired-end data, and the --poly_g_min_len option is used to specify the minimum length of poly-G tail to trim. The -h option is used to specify the name of the report file.
