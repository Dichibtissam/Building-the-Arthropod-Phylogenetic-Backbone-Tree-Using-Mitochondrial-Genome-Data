#!/bin/bash

# Activate the conda environment
source /opt/src/miniforge3/etc/profile.d/conda.sh
conda activate macse

# Define input directory containing alignments (.fas files)
ALIGN_DIR=~/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/alignment_Macse

# Define output directory for trimmed files
OUTPUT_DIR=~/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/Trim_alignment_macse
mkdir -p "$OUTPUT_DIR"

# Loop through all .fas files in the input directory
for FILE in "$ALIGN_DIR"/*.fas; do
    BASENAME=$(basename "$FILE" .fas)

    # Define output filenames
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_trimmed.fas"
    STATS_FILE="${OUTPUT_DIR}/${BASENAME}_trim_stats.csv"

    # Run trimming with MACSE
    macse -prog trimAlignment -align "$FILE" -respect_first_RF_ON \
          -min_percent_NT_at_ends 0.8 -half_window_size 3 \
          -out_trim_info "$STATS_FILE" \
          -out_NT "$OUTPUT_FILE"

    echo "✅ Trimmed: $FILE -> $OUTPUT_FILE"
done

