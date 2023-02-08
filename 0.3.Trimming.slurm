##################################################################################################################################
##### This is the command we used for filtering, here you can see the filter for trimming adapters and trimming of poliG > 2 #####
##################################################################################################################################

module load bioinfo-tools
module load fastp
i=/path/to/R1_001.fastq.gz
j=/path/to/R2_001.fastq.gz
h=/path/to/results/folder
g=MeDIP.R1.nopoliG.fastq.gz #name of the R1 output
f=MeDIP.R2.nopoliG.fastq.gz #name of the R2 output
w=MeDIP.nopoliG.html #name of the report

fastp --in1 $i --in2 $j --out1 $h/$g --out2 $h/$f -z 4 --detect_adapter_for_pe --poly_g_min_len 2 -h $h/$w
