
#!/bin/bash

################################# Loop the process for automated download, QC & Salmon Quantification ####################


################################# DEFINE DIRECTORY PATHS ###################################

data="../data"
process_data="../data/processed_reads"
acc_list="../data/SRR_Acc_List.txt"
ref_files="../ref_files"
fastqc="../fastqc"
fastp="../fastp_output"
salmon_output="../salmon_output"



##################################### LOOP THE PROCESS ######################################


########## Start from 3rd accession defined in txt file ##########################

tail -n +1 "$acc_list" | while read -r accession; do
   
    echo "Processing Accession: $accession"


    ################################# STEP-1: DOWNLOAD FASTQ FILES #######################################

    echo "Starting Download: $accession"

    fasterq-dump "$accession" --outdir "${data}" --progress 

    echo "Download Complete!"

    ################################### STEP-2: COMPRESS FILES TO SAVE DISK SPACE ###############################

    echo "Start Compressing: $accession"

    gzip ${data}/$accession*.fastq

    echo "Compression Complete!"

    ################################### STEP-3: QUALITY CHECK:FASTQC ##############################################

    echo "Start Fastqc: $accession"

    fastqc -o ${fastqc} ${data}/$accession*.fastq.gz

    echo "QC Complete!"

    ###################################### STEP-4: QUALITY CONTROL:FASTP ###############################################

    echo "Start fastp: $accession"

    fastp -i ${data}/${accession}_1.fastq.gz \
          -I ${data}/${accession}_2.fastq.gz \
          -o ${process_data}/${accession}_1.fastq.gz \
          -O ${process_data}/${accession}_2.fastq.gz \
          -h ${fastp}/${accession}_fastp.html \
          -j ${fastp}/${accession}_fastp.json

    echo "Fastp Complete!"

    ########################################### STEP-5: SALMON QUANTIFICATION #########################################

    echo "Salmon Quantification: $accession"

    salmon quant -i ${ref_files}/salmon_index -l A \
        -1 ${process_data}/${accession}_1.fastq.gz \
        -2 ${process_data}/${accession}_2.fastq.gz \
        -o ${salmon_output}/${accession}_trans_quant \
        --validateMappings 

    echo "Quantification Complete!"

    ######################################### STEP-6: CHECK MAPPING RATE ######################################

    cat ${salmon_output}/${accession}_trans_quant/logs/salmon_quant.log | grep "Mapping rate"

    ###################################### STEP7: RENAME QUANT FILE #####################################
    mkdir -p ${salmon_output}/quant_files/${accession}

    mv ${salmon_output}/${accession}_trans_quant/quant.sf \
       ${salmon_output}/quant_files/${accession}/quant.sf

    echo "Finished processing: $accession."

    ######################################## Continue with next accession ###########################################

read -p "Do you want to continue with the next accession? (y/n): " confirm < /dev/tty
    if [[ "$confirm" != "y" ]]; then
        echo "Stop Processing"
        break
    fi
done



##################################### END OF SCRIPT ##########################################################

