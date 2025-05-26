#!/bin/bash

# Usage:
# ./select_complete_representatives.sh tree_file n alignments_dir output_dir

# Prepare conda environment
source /opt/src/miniconda3/etc/profile.d/conda.sh
conda activate parnas_env

set -euo pipefail

TREE_FILE="$1"
N_REPRESENTATIVES="$2"
ALIGNMENTS_DIR="$3"
OUTPUT_DIR="$4"

TMPDIR=$(mktemp -d)
COMMON_IDS="$TMPDIR/common_ids.txt"

echo "🔍 Searching for IDs present in all gene files..."

# Extract accession IDs from each gene file (remove gene suffix after "_")
for file in "$ALIGNMENTS_DIR"/*.fas; do
    gene_name=$(basename "$file" .fas | sed 's/_cleaned$//')
    grep '^>' "$file" | sed 's/>//' | sed "s/_$gene_name\$//" > "$TMPDIR/$(basename "$file" .fas).ids"
done

# Get intersection of IDs that are present in at least 90% of the gene files
total_files=$(ls "$TMPDIR"/*.ids | wc -l)
threshold=$(echo "$total_files * 0.9" | bc | awk '{print int($1+0.9999)}')  # round up

ls "$TMPDIR"/*.ids | xargs cat | sort | uniq -c | awk -v threshold="$threshold" '$1 >= threshold {print $2}' > "$COMMON_IDS"

n_common_ids=$(wc -l < "$COMMON_IDS")
echo "✔️ Common IDs found: $n_common_ids"

if [ "$n_common_ids" -lt 2 ]; then
    echo "❌ ERROR: Not enough common taxa across gene alignments. Exiting."
    exit 1
fi

# Ensure requested number of representatives is less than number of taxa
if [ "$N_REPRESENTATIVES" -ge "$n_common_ids" ]; then
    echo "⚠️ WARNING: Requested number of representatives (n=$N_REPRESENTATIVES) >= number of available taxa ($n_common_ids)."
    echo "👉 Reducing n to $(($n_common_ids - 1))"
    N_REPRESENTATIVES=$(($n_common_ids - 1))
fi

# Build cleaned tree with only common IDs
echo "🌲 Pruning tree to keep only common taxa..."
pruned_tree="$TMPDIR/pruned_tree.nwk"
nw_prune "$TREE_FILE" -v -f "$COMMON_IDS" > "$pruned_tree"

# Run Parnas to select n representatives from the pruned tree
echo "🧬 Running Parnas to select $N_REPRESENTATIVES representatives..."
mkdir -p "$OUTPUT_DIR"
parnas_output="$OUTPUT_DIR/representatives.txt"
parnas -t "$pruned_tree" -n "$N_REPRESENTATIVES" > "$parnas_output"

echo "✅ Done. Selected representatives saved in: $parnas_output"

# Extract sequences of selected representatives from each gene fasta into cleaned fasta files
echo "📂 Extracting sequences of selected representatives for each gene..."

# Prepare pattern file: accession numbers with underscore (e.g., "^>OP156999_")
awk '{print "^>"$1"_"}' "$parnas_output" > "$TMPDIR/rep_patterns.txt"
head -n 5 "$ALIGNMENTS_DIR"/*.fas

# Extract sequences of selected representatives from each gene fasta into cleaned fasta files
echo "📂 Extracting sequences of selected representatives for each gene..."

for fasta in "$ALIGNMENTS_DIR"/*.fas; do
    gene_name=$(basename "$fasta" .fas | sed 's/_cleaned$//')
    output_fasta="$OUTPUT_DIR/${gene_name}_cleaned_representatives.fas"
    echo "🧬 Processing gene: $gene_name"

    awk -v ids="$parnas_output" -v gene="$gene_name" '
        BEGIN {
            while ((getline line < ids) > 0) {
                gsub(/\r/, "", line);  # remove carriage returns
                wanted[line"_"gene] = 1
            }
            RS = ">"; ORS = "";
        }
        NR > 1 {
            header = substr($0, 1, index($0, "\n")-1)
            seq = substr($0, index($0, "\n")+1)
            if (header in wanted) {
                print ">"$0
            }
        }
    ' "$fasta" > "$output_fasta"
done

echo "✅ All representative sequences extracted."

# Clean up
rm -r "$TMPDIR"

