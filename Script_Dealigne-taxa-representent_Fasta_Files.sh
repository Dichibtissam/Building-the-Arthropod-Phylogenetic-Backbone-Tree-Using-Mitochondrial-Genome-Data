#!/bin/bash

# Input directory containing the raw concatenated gene FASTA files
input_dir="/home/idich/internship/ibtissam/tree/AllGroups_GeneWise_RAW"

# Output directory for the de-aligned gene FASTA files
output_dir="/home/idich/internship/ibtissam/tree/AllGroups_GeneWise_DESALIGNED"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .fas files in the input directory
for file in "$input_dir"/*.fas; do
  echo "De-aligning $(basename "$file")"

  # Use seqkit to remove gaps and write to the output directory
  seqkit seq --remove-gaps "$file" > "$output_dir/$(basename "$file")"

  echo "  → Saved: $output_dir/$(basename "$file")"
done

echo "✔️ All gene files have been successfully de-aligned."

