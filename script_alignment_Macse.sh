#!/bin/bash

# Activate the conda environment
source /opt/src/miniforge3/etc/profile.d/conda.sh
conda activate macse

# Define input and output directories
INPUT_DIR="/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/"
OUTPUT_DIR="/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/alignment_Macse"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# File to store the list of jobs
JOB_LIST="list-job"
rm -f $JOB_LIST  # Remove the old job list if it exists

# Generate the list of alignment commands to run in parallel
echo "ls ${INPUT_DIR}/*.fas"
for SEQ in $(ls ${INPUT_DIR}/*.fas); do
    echo $SEQ
    GENE_NAME=$(basename "$SEQ" .fas)
    OUTPUT_FILE="${OUTPUT_DIR}/${GENE_NAME}.fas"

    # Add alignment command to the job list
    echo "macse -prog alignSequences -seq $SEQ -gc_def 5 -out_NT $OUTPUT_FILE" >> $JOB_LIST
done

# Run the commands in parallel (10 simultaneous tasks)
cat $JOB_LIST | parallel -j 10 --eta
echo "✅ Alignment completed for all genes in parallel."



