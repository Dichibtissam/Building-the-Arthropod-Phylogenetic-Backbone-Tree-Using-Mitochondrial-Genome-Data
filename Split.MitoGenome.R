#Function developed by David Eme
#' @title function to Split the wholemitogenomes into x different fasta alignment belonging to the x gene region of interests.  

#' @description This function exports a separate alignement in fasta format for each gene region of interest from a single long DNA sequence (e.g. a mitogenome).

#' @details This function check which accession number have all the gene of interest and export in fasta only the sequences of those accession numbers.

#' @param input.fasta an "alignment" object imported using the read.alignment function from the sequinr package.

#' @param Gene.to.get a vector listing all the gene name that should be exported.

#' @param Annotations a data frame containing the accession nimber the gene name of interest and the startig and ending position of thebgene on the fasta alignment.

#' @param Annotation.col the number of the colomn providing the gene name on the Annotations data frame.

#' @param Access.NB.co the number of the colomn providing the accession number on the Annotations data frame.

#' @param Starting.col the number of the colomn providing the starting position of the sequence on the Annotations data frame.

#' @param Ending.col the number of the colomn providing the ending position of the sequence on the Annotations data frame.

#' @param output.folder the path toward the output folder storing the fasta alignment fo the gene of interests.

#' @param return The function return the fasta alignments of the differents gene of interests in the output folder. It also returns into the R environment the annotation table including only the Accession numbers with all the genes of interest. 

require(seqinr)


Split.MitoGenome = function(input.fasta = NULL, Gene.to.get = NULL, Annotations = NULL, Annotation.col = NULL, Access.NB.col = NULL, Starting.col = NULL, Ending.col = NULL, output.folder = NULL){

### Assess which accession numbers have a complete set of gene.to.get.

Acces.NB = unique(Annotations[,Access.NB.col])

Complete.Mito.Acces.List = unlist(lapply(1:length(Acces.NB), function(x){
a = Annotations[which(Annotations[,Access.NB.col] == Acces.NB[x]),]
if(length(na.omit(match(Gene.to.get, a[,Annotation.col]))) == 13){
Acces.NB[x]
}
}))

### Get the subset of Accession with the complete set of 13 genes.
# Annot.sub = Annotations[match(Complete.Mito.Acces.List, Annotations[,Access.NB.col]),]
Annot.sub = Annotations[Annotations[,Access.NB.col] %in% Complete.Mito.Acces.List,]


Annot.sub$Access.Nb = unlist(lapply(strsplit(Annot.sub$Access.Nb, ".", fixed = T), function(j){
j[1]}))

# Remove the version number of the accession number.
b = unlist(lapply(strsplit(Complete.Mito.Acces.List, ".", fixed = T), function(j){
j[1]}))

# Export the different alignment
lapply(1:length(Gene.to.get), function(x){
Align.file = paste(output.folder, "/", Gene.to.get[x], ".fas", sep = "")
cat(file = Align.file)
lapply(1:length(b), function(i){

Annot.sub.Gene = Annot.sub[which(Annot.sub[,Annotation.col] == Gene.to.get[x]),]

starting = as.numeric(as.character(Annot.sub.Gene[which(Annot.sub.Gene$Access.Nb == b[i]), Starting.col]))
ending = as.numeric(as.character(Annot.sub.Gene[which(Annot.sub.Gene$Access.Nb == b[i]), Ending.col]))

# Check if the starting position is alwaysthe first or if the two are inverted.
if(ending < starting){
sequence.extract = substr(input.fasta$seq[which(input.fasta$nam == b[i])], ending, starting)
} else {
sequence.extract = substr(input.fasta$seq[which(input.fasta$nam == b[i])],  starting, ending)
}

cat(file = Align.file, ">",b[i], "_", Gene.to.get[x], "\n", sequence.extract , "\n", sep = "", append = T)
})
})


# Corret the Annot.sub
i = 1
for( i  in 1:dim(Annot.sub)[1]){
starting = as.numeric(as.character(Annot.sub[i,Starting.col]))
ending = as.numeric(as.character(Annot.sub[i, Ending.col]))

if(as.numeric(as.character(ending)) < as.numeric(as.character(starting))){
Annot.sub[i,Starting.col] = ending
Annot.sub[i,Ending.col] = starting
}
}

return(Annot.sub)
}

# For Debug
# rm(input.fasta, a, b, Gene.to.get, Annotations, Annotation.col, Access.NB.col, Starting.col, Ending.col, output.folder, Align.file, sequence.extract, starting, ending, Annot.sub, Complete.Mito.Acces.List, Acces.NB)

