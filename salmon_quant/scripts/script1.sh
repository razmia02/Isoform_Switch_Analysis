#!/bin/bash

##################### Script to download fastq files, initial QC, Salmon Index & Salmon Quantification ###########################

#_____________________________________________________________________________________________________________________________________

#################################### SETUP DIRECTORIES ######################################################

mkdir data fastqc salmon_output scripts ref_files


#################################### DEFINE DIRECTORY PATHS ###########################################

data=~/salmon_quant/data
process_data=~/salmon_quant/data/processed_reads
acc_list=/mnt/d/Projects/Isoform_Switching/SRR_Acc_List.txt
ref_files=~/salmon_quant/ref_files
fastqc=~/salmon_quant/fastqc
fastp=~/salmon_quant/fastp_output
salmon_output=~/salmon_quant/salmon_output

#___________________________________________________________________________________________________________________

#################################### SETUP PACKAGES & LIBRARIES  ####################################################

######################## SRA-TOOLKIT ##############################

sudo apt install sra-toolkit ####### Install Sra-toolkit

fastq-dump --version #### Check if sra-toolkit is installed successfully

###################### SALMON ##############################

### Add the necessary channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

############### Install Salmon ##########################

conda create -n salmon -c bioconda -c conda-forge salmon

conda activate salmon

salmon --version

#_______________________________________________________________________________________________________________________________

################################## DOWNLOAD REFERENCE FILES  ######################################################

############################# Reference Transcripts from Gencode ###########################

wget -P ${ref_files} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz 

################################# GTF file for annotation ##################################################

wget -P ${ref_files} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz


######################### Genome Fasta for decoy aware quantification #########################################

wget -P ${ref_files} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz


#__________________________________________________________________________________________________________________________

########################################### STEP-1: DOWNLOAD FASTQ FILES  ############################################

#################### NCBI GEO DATASETID = GSE228870 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228870)

################### 7 control & 7 tumor samples have been used for this project ##############################

############## The list of accessions used are in txt file (SRR_Acc_List.txt) ##########################

################# Download the file ##########################

echo "Starting Download"

fasterq-dump SRR24053380 --outdir ${data} --progress

echo "Download Complete!"

########## Compress the files to save disk space #############

echo "Compressing Files"

gzip ${data}/SRR24053380*.fastq

echo "Compressed"

#___________________________________________________________________________________________________________________________

############################################ STEP-2: QUALITY CHECK ##########################################################

######################### Run FASTQC on the files ######################################

echo "Start QC"

fastqc -o ${fastqc} ${data}/SRR24053380*.fastq.gz

echo "QC Complete!"

#_____________________________________________________________________________________________________________________________________

####################################### STEP-3: QUALITY CONTROL ####################################################################

########################## FASTP ################################

echo "Start fastp"

fastp -i ${data}/SRR24053380_1.fastq.gz -I ${data}/SRR24053380_2.fastq.gz -o ${process_data}/SRR24053380_1.fastq.gz \
        -O ${process_data}/SRR24053380_2.fastq.gz -h ${fastp}/SRR24053380_fastp.html

echo "Fastp Complete!"


########### Run Fastqc again to check the quality of processed reads ##################

#______________________________________________________________________________________________________________________________
 
################################################ STEP-4: SALMON INDEX ############################################################

echo "Starting Index"

salmon index -t ${ref_files}/gencode.v44.transcripts.fa.gz -i ${ref_files}/salmon_index -k 31 --gencode

echo "Index Complete!"

#__________________________________________________________________________________________________________________________________

####################################### STEP-5: SALMON QUANTIFICATION #################################################################

echo "Salmon Quantification"

salmon quant -i ${ref_files}/salmon_index -l A \
        -1 ${process_data}/SRR24053380_1.fastq.gz -2 ${process_data}/SRR24053380_2.fastq.gz  \
        -o ${salmon_output}/SRR24053380_trans_quant --validateMappings

echo "Quantification Complete!"

cat ${salmon_output}/SRR24053380_trans_quant/logs/salmon_quant.log | grep "Mapping rate" ######### Check the mapping rate

mv ${salmon_output}/SRR24053380_trans_quant/quant.sf ${salmon_output}/SRR24053380_trans_quant/SRR24053380_quant.sf ######### Rename quant.sf file

echo "Process Completed!"


########################### Check if transcripts match between annotation, transcripts and quant.sf files #################

# Get the first transcript from the FASTA
zgrep ">" ${ref_files}/gencode.v44.transcripts.fa.gz | head -n 1 | cut -d'|' -f1

# Get the same transcript from the GTF
zgrep "ENST00000456328" ${ref_files}/gencode.v44.annotation.gtf.gz | head -n 1 | grep -o "transcript_id \"[^\"]*\""

grep "ENST00000456328" ${salmon_output}/SRR24053380_trans_quant/SRR24053380_quant.sf

#################################### quant.sf files will be used for further analysis #############################################

