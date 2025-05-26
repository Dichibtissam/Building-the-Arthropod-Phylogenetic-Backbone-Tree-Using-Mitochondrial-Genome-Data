library(ape)

# Define input directories for each group
input_dirs <- c(
  "~/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Anostraca/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Orthoptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Ixodida/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Decapoda/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Poduromorpha/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Spirostreptida/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Mesostigmata/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Notostraca/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Araneae/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Blattodea/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Stomatopoda/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Plecoptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Balanomorpha/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Scorpiones/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Trombidiformes/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Mantodea/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Entomobryomorpha/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Megaloptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Neuroptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Ephemeroptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Sarcoptiformes/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Odonata/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Amphipoda/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Mecoptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Phasmatodea/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Raphidioptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Phthiraptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Thysanoptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Siphonaptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Trichoptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Anaspidacea/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Hemiptera/Results/Separate.Align/Tree/File_backup",
  "~/internship/ibtissam/alignement/Diptera/Results/Separate.Align/Tree/File_backup"
)

# Output directory for merged gene-wise FASTA files
output_dir <- "~/internship/ibtissam/tree/AllGroups_GeneWise_RAW"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# List of mitochondrial genes of interest
genes <- c("atp6", "atp8", "co1", "co2", "co3", "cytb",
           "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6")

# Function to read a FASTA file and return named sequences
read_fasta_safely <- function(file_path) {
  if (file.exists(file_path)) {
    readLines(file_path)
  } else {
    NULL
  }
}

# Loop through each gene
for (gene in genes) {
  cat("Processing gene:", gene, "\n")
  all_lines <- character()

  for (dir in input_dirs) {
    file_path <- file.path(dir, paste0(gene, "_cleaned_representatives.fas"))
    fasta_lines <- read_fasta_safely(file_path)

    if (!is.null(fasta_lines)) {
      all_lines <- c(all_lines, fasta_lines)
    } else {
      cat("  ❌ Missing file in:", dir, "\n")
    }
  }

  # Write combined FASTA file
  if (length(all_lines) > 0) {
    output_file <- file.path(output_dir, paste0(gene, "_ALL_raw.fas"))
    writeLines(all_lines, output_file)
    cat("  ✅ Written to:", output_file, "\n")
  } else {
    cat("  ⚠️ No sequences found for gene:", gene, "\n")
  }
}
