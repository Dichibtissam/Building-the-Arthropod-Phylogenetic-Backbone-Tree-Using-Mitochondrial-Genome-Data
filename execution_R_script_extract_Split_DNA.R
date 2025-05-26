# lines before 120...

# # Load the required libraries
library(rentrez)  # To interact with the GenBank database                                                       
library(XML)      # To process and extract XML information from GenBank records

# Load the data from the CSV file containing information about the "in" and "out" groups
file_path <- "outgroup_selection_with_accessions_and_taxids.csv"


group_data <- read.csv(file_path)  # Read the CSV file into a variable 'group_data'

# Filter the row containing "Diplostraca" in the "order" column
Diplostraca_data <- group_data[grepl("Diplostraca", group_data$order, ignore.case = TRUE), ]

# Check that we do have a row with "Diplostraca"
print(Diplostraca_data)

# Extract accession numbers from the "ingroup_accessions" and "outgroup_accessions" columns
ingroup_accessions <- Diplostraca_data$ingroup_accessions
outgroup_accessions <- Diplostraca_data$outgroup_accessions

# Merge the two accession lists into a single list
Test.Accession.Diplostraca <- unique(c(ingroup_accessions, outgroup_accessions))

# Split the elements if they contain commas (in case some entries are comma-separated strings)
Test.Accession.Diplostraca <- unlist(strsplit(Test.Accession.Diplostraca, ","))

# Remove unnecessary spaces around the accession names
Test.Accession.Diplostraca <- trimws(Test.Accession.Diplostraca)  # Clean up any spaces

# Remove duplicates again in case the splitting generated redundancies
Test.Accession.Diplostraca <- unique(Test.Accession.Diplostraca)

# Display the accession numbers
print(Test.Accession.Diplostraca)

# Check if there are any NA (missing) values in the vector Test.Accession.Diplostraca
which(is.na(Test.Accession.Diplostraca)) # no NA reported.

# Search for strings that contain "NA" as text (not actual missing values)
Test.Accession.Diplostraca[grep("NA", Test.Accession.Diplostraca)] # no NA reported

#######################
### Extract data Diplostraca

# Load the function to extract the mitogenome in a fasta file, as well as the meta data and the sequence annotations                 
source("Scripts/Diplostraca/Get.Seq.Metadata.NCBI.R")

####
# Record the start time of the process to measure execution time
cpu0 = Sys.time()

# Retrieve sequence metadata from the NCBI database
# using the accession numbers stored in 'Test.Accession.Diplostraca'
# and store the results in the object 'Data.Diplostraca'
# The 'output.Path' parameter specifies where to save the retrieved files
# and names the output file "Test.Extraction.Diplostraca.And.outgroup"
Data.Diplostraca = Get.Seq.Metadata.NCBI(Accession.NB = Test.Accession.Diplostraca, output.Path = "/home/idich/internship/ibtissam/alignement/Diplostraca", output.Name = "Test.Extraction.Diplostraca.And.outgroup")

# Record the end time of the process
cpu1 = Sys.time()

# Calculate and display the total execution time for data retrieval
cpu1 - cpu0 #13.98732 mins

# Display the column (or component) names of the 'Data.Diplostraca' object,
# which is a list or data frame containing the extracted metadata
names(Data.Diplostraca)

## Homogeneization of the gene name

# Convert gene names to lowercase to standardize them
Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = tolower(Data.Diplostraca$Annotation.Sequence.DF$gene_name)

# Replace "ox" with "o" in the corrected gene names
Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("ox", "o", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("iii", "3", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("ii", "2", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("i", "1", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("nadh", "nd", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("nad", "nd", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub(" ", "", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("atpase", "atp", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("cob", "cytb", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("co21", "co3", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("ndh-u", "nd", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("atp6gene", "atp6", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("atp8gene", "atp8", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("mt-", "", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("cyb", "cytb", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("cb", "cytb", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("ctyb", "cytb", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("nd-", "nd", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

# Display how many genes are associated with each corrected gene name
table(Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)


##########################################################################################
#### Export the 13 genes of interests in 13 fsata file : emprical example with the Diplostraca.

# lines after line 120 are commented out
require(seqinr)

# Load the Mitogenome exported by the function "Get.Seq.Metadata.NCBI.R" into R in as an "alignment" object
Mitogenomes.fasta = seqinr::read.alignment("/home/idich/internship/ibtissam/alignement/Diplostraca/Test.Extraction.Diplostraca.And.outgroup.Align.fas",format = "fasta")

# source the function tosplit the mitogeoes into 13 fasta alignment
source("Scripts/Diplostraca/Split.MitoGenome.R")

# Create a "Results" directory in the specified path to store the alignment results
dir.create("/home/idich/internship/ibtissam/alignement/Diplostraca/Results")

# Create a subdirectory "Separate.Align" to organize the alignment results of the separated genes
dir.create("/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align")

# Run the function, export the 13 gene alignments in fasta format.
# see directly n the file of the function to get the information about the different arguments.ø
Annot.DF.Diplostraca.test = Split.MitoGenome(input.fasta = Mitogenomes.fasta, Gene.to.get = c("co1", "co2", "co3", "cytb", "atp6", "atp8", "nd1", "nd2", "nd3","nd4", "nd4l", "nd5", "nd6"), Annotations = Data.Diplostraca$Annotation.Sequence.DF, Annotation.col = 6, Access.NB.col = 1, Starting.col = 3, Ending.col = 4, output.folder = "/home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align")

# Exit the R session
quit()


# Change directory to move to the "Separate.Align" folder
cd /home/idich/internship/ibtissam/alignement/Diplostraca/Results/Separate.Align

# Check the presence of outgroup accession numbers in "atp6.fas"
grep -E "KU058637|MN356346|MT017889|MW393772|JQ071617" atp6.fas

# Count the number of sequences in a FASTA file (where each sequence starts with a ">")
grep -c ">" atp6.fas







