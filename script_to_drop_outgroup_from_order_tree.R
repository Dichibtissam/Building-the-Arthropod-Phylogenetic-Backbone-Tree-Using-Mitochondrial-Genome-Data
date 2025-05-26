#!/usr/bin/env Rscript

library(ape)
library(dplyr)
library(stringr)

# Load CSV metadata with accessions and orders
meta <- read.csv("/home/idich/internship/ibtissam/alignement/outgroup_selection_with_accessions_and_taxids.csv", stringsAsFactors = FALSE)

# Define the target order to filter
target_order <- "Diplostraca"

# Extract raw accession numbers for that order
accessions_raw <- meta %>%
  filter(order == target_order) %>%
  pull(outgroup_accessions)

# Remove suffixes like ".1", ".2" from accession numbers
# Split by comma, remove whitespace, flatten, and clean suffixes
accessions_clean <- accessions_raw %>%
  str_split(",\\s*") %>%         # split by comma and optional space
  unlist() %>%
  str_replace("\\.\\d+$", "")    # remove suffixes like ".1"

# Print the cleaned list of outgroup accessions
cat("Outgroup accessions to remove:\n")
print(accessions_clean)


# Path to the original tree file
tree_path <- "/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/Tree/File_backup/Partitions_Concat_3codons.nex.best_scheme.nex.contree"

# Read the tree
tree <- read.tree(tree_path)


# Check which accessions exist in the tree
accessions_in_tree <- accessions_clean[accessions_clean %in% tree$tip.label]
accessions_not_found <- setdiff(accessions_clean, tree$tip.label)

# Message about found and not found accessions
cat("\nAccessions found in tree:\n")
print(accessions_in_tree)

cat("\nAccessions NOT found in tree:\n")
print(accessions_not_found)


cat(sprintf("\nTotal to remove: %d | Found in tree: %d | Not found: %d\n",
            length(accessions_clean),
            length(accessions_in_tree),
            length(accessions_not_found)))


# Drop the tips that correspond to the cleaned accession numbers
tree2 <- drop.tip(tree, accessions_in_tree)


# Construct output path in the same folder as input tree
output_path <- file.path(dirname(tree_path), "Partitions_Concat_3codons_no_outgroup.nwk")

# Save the pruned tree in Newick format to the output path
write.tree(tree2, file = output_path)


cat("\nFiltered tree saved to:\n", output_path, "\n")
