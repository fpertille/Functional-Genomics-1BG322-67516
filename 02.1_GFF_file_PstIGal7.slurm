#!/bin/bash -l
#SBATCH -A uppmax2023-2-10
#SBATCH -M snowy
#SBATCH -p core -n 1
#SBATCH --mem=64gb
#SBATCH -t 24:00:00
#SBATCH -J trimm
#SBATCH --error /error_path/00.PreparingInSilicoFile.err
#SBATCH --output /output_path/00.PreparingInSilicoFile.out

############################################################################################
#### the output from 02.in_silico_digestion_PstIGal7.R should be uploaded to the server ####
############################################################################################

#path to the folder
i=/path/to/folder
#output from 02.in_silico_digestion_PstIGal7.R
h=in_silico_windows_galgal7.gff
#code for the reference genome in use
j=galgal7

# extracting the header of the gff file
head -n 3 $i/$h > $i/header.txt

# get the gff file without the header
sed '1,3d' $i/$h > $i/raw_windows.txt

#Creating a file with the ID name for as many times as cut you have in the ggf file | COLUMN 9
START=1
END=$(wc -l < $i/raw_windows.txt)
echo $END
for (( c=$START; c<=$END; c++ ))
do
    echo "ID=$c;gene_id=$c;"
done > $i/columns9.txt

paste $i/raw_windows.txt $i/columns9.txt >> $i/header.txt
mv $i/header.txt $i/in_silico_windows.curated.$j.gff
rm $i/raw_windows.txt
#rm $i/header.txt
rm $i/columns9.txt
