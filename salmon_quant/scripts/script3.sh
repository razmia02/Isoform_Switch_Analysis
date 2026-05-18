#!/bin/bash

###################################### SCRIPT TO RUN PFAM Locally #####################################################

############ Specify Directories ####################


down_analysis=~/salmon_quant/downstream_analysis
fasta_file_nt=/mnt/d/Projects/Isoform_Switching/isoformswitchanalyzer/switchanalyzer_output/isoformSwitchAnalyzeR_isoform_nt.fasta
fasta_file_aa=~/salmon_quant/downstream_analysis/isoformSwitchAnalyzeR_isoform_AA.fasta
pfam_output=~/salmon_quant/downstream_analysis

######################### Install pfam scan ########################

conda install -c bioconda hmmer -y

############### Check Version #################

hmmscan -h 

###################### Download required files ###############################

wget -O ${down_analysis}/Pfam-A.hmm.dat.gz http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz

wget -O ${down_analysis}/Pfam-A.hmm.gz http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

################## Unpack the downloaded files ############################3


gunzip -c ${down_analysis}/Pfam-A.hmm.dat.gz > ${down_analysis}/Pfam-A.hmm.dat

gunzip -c ${down_analysis}/Pfam-A.hmm.gz > ${down_analysis}/Pfam-A.hmm

rm ${down_analysis}/Pfam-A.hmm.gz ${down_analysis}/Pfam-A.hmm.dat.gz


####################### Prepare pfam database for HMMER by creatin binary files ###########################

hmmpress ${down_analysis}/Pfam-A.hmm


######################### Run hmmscan  ##############################

hmmscan --domtblout ${down_analysis}/pfam_results.txt ${down_analysis}/Pfam-A.hmm ${fasta_file_aa}


########################## Copy the pfam_results.txt to your RStudio Project directory ############################################