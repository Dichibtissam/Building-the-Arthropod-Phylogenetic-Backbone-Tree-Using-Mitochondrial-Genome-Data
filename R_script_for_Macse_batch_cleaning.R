## Loading the necessary libraries
require(regPhylo)  # Pour le nettoyage et l'analyse des séquences alignées
require(ape)       # Pour la manipulation des séquences ADN

## Loading the necessary libraries
alignment_dir <- "/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/alignment_Macse/"
alignment_files <- list.files(alignment_dir, pattern = "\\.fas$", full.names = TRUE)

## Create a directory for output files (if not existing)
output_dir <- "/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align/cleaned_Macse_alignment"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Create a report file
report_file <- file.path(output_dir, "cleaning_report.txt")
cat("Rapport de nettoyage des alignements\n", file = report_file)

## Loop to process each alignment file
for (file in alignment_files) {
  cat("\nTraitement du fichier :", file, "\n")
  cat("\nTraitement du fichier :", file, "\n", file = report_file, append = TRUE)
  
  ## Load the alignment in FASTA format with error handling
  alignment_data <- tryCatch(
    ape::read.dna(file, format = "fasta"),
    error = function(e) {
      cat("Erreur lors de la lecture du fichier :", file, "\n")
      cat("Message d'erreur:", e$message, "\n")
      cat("Erreur lors de la lecture du fichier :", file, ":", e$message, "\n", file = report_file, append = TRUE)
      return(NULL)
    }
  )
  
  if (is.null(alignment_data)) {
  cat("Le fichier", file, "n'a pas pu être chargé.\n")
  cat("Le fichier", file, "n'a pas pu être chargé.\n", file = report_file, append = TRUE)
  next
}

# Check if the alignment data is empty
if (length(alignment_data) == 0) {
  cat("Le fichier", file, "est vide ou ne contient pas de séquences valides.\n")
  cat("Le fichier", file, "est vide ou ne contient pas de séquences valides.\n", file = report_file, append = TRUE)
  next
}

  ## Remove problematic sequences (insertions, short sequences, ambiguities, etc.) with error handling
  cleaned_alignment <- tryCatch(
    rm.del.gap(
      input = alignment_data,
      Remove.del.gaps = TRUE,
      Remove.Seq.Indels = TRUE,
      Min.Perc.Seq.With.Info = 0.05,
      Nb.Seq.Indels = 1,
      Remove.Short.Seq = FALSE,
      Max.Per.missing = 85,
      Remove.Seq.TooManyAmbig = FALSE,
      Percent.Ambig = 30,
      Remove.Loc.Low.Freq = FALSE,
      Minimal.Locus.Freq = 99
    ),
    error = function(e) {
      cat("Erreur lors du nettoyage du fichier :", file, "\n")
      cat("Message d'erreur:", e$message, "\n")
      cat("Erreur lors du nettoyage du fichier :", file, ":", e$message, "\n", file = report_file, append = TRUE)
      return(NULL)
    }
  )
  
  if (is.null(cleaned_alignment)) next
  
  ## Check and log removed sequences
  removed_sequences <- if (!is.null(cleaned_alignment$Sequences.Removed) && length(cleaned_alignment$Sequences.Removed) > 0) {
    paste(unlist(cleaned_alignment$Sequences.Removed), collapse = ", ")
  } else {
    "Aucune"
  }
  cat("Séquences supprimées :", removed_sequences, "\n")
  cat("Séquences supprimées :", removed_sequences, "\n", file = report_file, append = TRUE)
  
  ## Detect stop codons with error handling
  stop_codon_check <- tryCatch(
    Detect.stop.codon(
      input = cleaned_alignment$Alignment,
      genet.code = 5
    ),
    error = function(e) {
      cat("Erreur lors de la détection des codons stop :", file, "\n")
      cat("Erreur lors de la détection des codons stop :", file, "\n", file = report_file, append = TRUE)
      return(NULL)
    }
  )
  
  detected_stop_codons <- if (!is.null(stop_codon_check) && length(stop_codon_check) > 0) {
    paste(unlist(stop_codon_check), collapse = ", ")
  } else {
    "Aucun"
  }
  cat("Codons stop détectés :", detected_stop_codons, "\n")
  cat("Codons stop détectés :", detected_stop_codons, "\n", file = report_file, append = TRUE)
  
  ## Save the cleaned file to the output directory
  output_file <- file.path(output_dir, gsub("\\.fas", "_cleaned.fas", basename(file)))
  ape::write.dna(cleaned_alignment$Alignment, file = output_file, format = "fasta")
  cat("Fichier nettoyé sauvegardé sous :", output_file, "\n")
  cat("Fichier nettoyé sauvegardé sous :", output_file, "\n", file = report_file, append = TRUE)
}

## End of script messages
cat("✔️ All files processed and saved in:", output_dir, "\n")
cat("✔️ Report on deleted sequences saved in:", report_file, "\n")


## End of script

