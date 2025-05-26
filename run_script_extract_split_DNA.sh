# Read the first 120 lines of the script

script_lines <- readLines("/home/idich/internship/ibtissam/alignement/Scripts/Diplostraca/execution_R_script_extract_Split_DNA.R", n = 120)

# Parse and execute these lines
eval(parse(text = script_lines))




# Read all lines of the script
script_lines <- readLines("/home/idich/internship/ibtissam/alignement/Scripts/Diplostraca/execution_R_script_extract_Split_DNA.R")

# Subset the lines from 121 to 145
script_lines_121_to_145 <- script_lines[121:145]


# Parse and execute these lines

for (line in script_lines_121_to_145) {
  tryCatch(
    eval(parse(text = line)),
    error = function(e) message("Erreur ignorée : ", e$message)
  )
}

# Exit the current R session
quit()

# Change directory to move to the "Separate.Align" folder
cd /home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align

# Check the presence of outgroup accession numbers in "atp6.fas"
grep -E "AB084514|OL757477|MK579381|MF496656|AY639934" atp6.fas

# Count the number of sequences in a FASTA file (where each sequence starts with a ">")
grep -c ">" atp6.fas


