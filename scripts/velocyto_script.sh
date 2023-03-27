#!/bin/bash

#------------------------------------------------------------------------------------------------------------
# job: Run velocyto CLI to get the spliced and unspliced .loom files! 
# date: 30/05/2021
# author: António Sousa

echo "---"
echo ""
echo "Job: Run velocyto CLI to get the spliced and unspliced .loom files!" 
echo "Date: `date +%d/%m/%Y`"
echo "Author: António Sousa - UBI-IGC"
echo "Script run at: `hostname` (IP: `hostname -I`)"
echo ""
echo "---"
echo ""
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
# abort if some command returns a non-zero value
set -e
set -o pipefail
#------------------------------------------------------------------------------------------------------------
## Variables 

# PATH 
AFS_PATH=/afs/igc.gulbenkian.pt/folders/UBI/PROJECTS/UBI-2021/ongoing/2012_miguel_elisa

# GTF and bam files for Cre3 and Lox2 samples
SAMPLES='Cre3 Lox2'
bam_cre3=${AFS_PATH}/data/ftp01.igc.gulbenkian.pt/Cre3_count_full/outs/possorted_genome_bam.bam
bam_lox2=${AFS_PATH}/data/ftp01.igc.gulbenkian.pt/Lox2_count_full/outs/possorted_genome_bam.bam
GTF=${AFS_PATH}/data/ftp01.igc.gulbenkian.pt/Pchabaudi2/genes/genes.gtf
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
## check input
echo ""
echo "---"
echo ""

# check if velocyto is installed
# activate conda first
source /home/agsousa/miniconda3/etc/profile.d/conda.sh
conda activate velocyto
if [[ -x "$(command -v velocyto)" ]];
then
        echo "";
        echo "velocyto is installed!";
        echo "$(velocyto --version)"
        echo "";
else
        echo "";
        echo "ABORT!";
        echo "velocyto is not installed or in your PATH!";
        exit 0
        echo ""
fi;

# check if GNU parallel is installed 
if [[ -x "$(command -v parallel)" ]];
then
        echo "";
        echo "GNU parallel is installed!";
        echo "$(parallel -V)"
        echo "";
else
        echo "";
        echo "ABORT!";
        echo "GNU parallel is not installed or in your PATH!";
        exit 0
        echo ""
fi;

# check folders
for folder in $AFS_PATH; 
do
        if [ ! -d $folder ]; 
        then 
                echo "Folder $folder does not exist! Exiting..."
                exit 0
        fi;
done;  

# check files
for file in $bam_cre3 $bam_lox2 $GTF;
do      
        if [ ! -f $file ]; 
        then
                echo "File $file does not exist! Exiting..."
                exit 0
        fi;     
done; 

echo ""
echo "---"
echo ""
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
## Run STAR: index and map
echo ""
echo "---"
echo ""
echo "Running velocyto CLI in order to get spliced/unspliced matrices..."
parallel --dry-run -v \
    velocyto run10x -@ 8 --samtools-memory 5000 \
    ${AFS_PATH}/data/ftp01.igc.gulbenkian.pt/{}_count_full \
    $GTF ::: $SAMPLES &> velocyto_script.log
SECONDS=0;
parallel -v \
    velocyto run10x -@ 8 --samtools-memory 5000 \
    ${AFS_PATH}/data/ftp01.igc.gulbenkian.pt/{}_count_full \
    $GTF ::: $SAMPLES >> velocyto_script.log 2>&1
conda deactivate
echo ""
echo "Finished velocyto!"
echo ""
echo "velocyto took $SECONDS secs to map reads!"
echo ""
echo "---"
echo ""
#------------------------------------------------------------------------------------------------------------
