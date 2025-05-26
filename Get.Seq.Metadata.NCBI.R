#' @title Extract DNA sequences,metadata and annotation from GenBank using accession number


#' @description This function extracts all the sequences, accession numbers, and related
#' information including annotation of the gene regions (CDS) from an GenBank accession number.
#' This function is dedicated to extract the whole mitochondrial genome and the CDS of the 13 genes.

#' @return This function returns a fasta alignment associated with the accession number, and a rds list of 4 tables, the first table reports the meta data of the DNA sequences, the second table reports the annotation of the sequences, the third table reports the accession number without a mitogenomes associated, and the fourth table reports the accession number of the mitochondiral genome with less than 13 mitochondrial genes. These four tables are also exported into the R environment.

#' @param Accession.NB a vector of the accession numbers

#' @param output.Path the name of the path toward the output folder.

#' @param output.Name name of the output files, the alignment uses the suffix ".Align.fas", and the rds files including the list of 4 tables uses the suffix ".MetaData.Annotations.rds".

require(rentrez)


Get.Seq.Metadata.NCBI = function(Accession.NB = NULL, output.Path = NULL, output.Name = NULL){

Align.name = paste(output.Path, "/",  output.Name, ".Align.fas", sep = "")                     
cat(file = Align.name)  


Meta.data.sequences = matrix(nrow  = length(Accession.NB), ncol = 14)
Annotation.Sequence.list = list()

x=1
for(x in 1:length(Accession.NB)){

oo = tryCatch(rentrez::entrez_search(db="nucleotide", 
				term = paste(Accession.NB[x], "[ACCN]", sep = ""),
                               use_history=TRUE), error=function(e) "Error")
                               
if(oo[1] == "Error"){

Meta.data.sequences[x,] = c(Accession.NB[x], rep(NA, 9), Accession.NB[x], NA, NA, NA)
Annot.seq.DF = as.data.frame(matrix(ncol = 4, nrow = 1))

} else {
                                  
oo.seqs = tryCatch(rentrez::entrez_fetch(db='nucleotide',
                                                web_history = oo$web_history, 
                                                rettype = 'gb', retmode = 'xml',
                                                parsed = TRUE), error=function(e) "Error")
                                                
}


if(oo.seqs[1] == "Error"){

Meta.data.sequences[x,] = c(Accession.NB[x], rep(NA, 9), Accession.NB[x], NA, NA, NA)
Annot.seq.DF = as.data.frame(matrix(ncol = 4, nrow = 1))


} else {
        
                                      
oo.seqs.xml = XML::xmlToList(oo.seqs)

if(class(oo.seqs.xml[[1]]) == "list"){

Meta.data.sequences[x,] = as.vector(unlist(oo.seqs.xml[[1]][match(c("GBSeq_locus", "GBSeq_length", "GBSeq_strandedness", "GBSeq_moltype", "GBSeq_topology", "GBSeq_division", "GBSeq_update-date", "GBSeq_create-date", "GBSeq_definition", "GBSeq_primary-accession", "GBSeq_accession-version", "GBSeq_source", "GBSeq_organism", "GBSeq_taxonomy"), names(oo.seqs.xml[[1]]))]))


cat(file = Align.name, ">",  oo.seqs.xml[[1]]$GBSeq_locus, "\n", oo.seqs.xml[[1]]$GBSeq_sequence, "\n", sep = "", append = T)


# Get the information about the 13 genes of the mitogenome, with the starting and ending position of each gene


feature.pos = match("GBSeq_feature-table", names(oo.seqs.xml[[1]]))
i = 1
Annot.seq.DF = as.data.frame(do.call(rbind, lapply(1:length(oo.seqs.xml[[1]][[feature.pos]]), function(i) {
#i=1
#for(i in 1:length(oo.seqs.xml[[1]][[feature.pos]])){
if(oo.seqs.xml[[1]][[feature.pos]][[i]][[1]] == "gene"){
Type = "gene"

interval.pos = match("GBFeature_intervals", names(oo.seqs.xml[[1]][[feature.pos]][[i]]))
starting = oo.seqs.xml[[1]][[feature.pos]][[i]][[interval.pos]][[1]]$GBInterval_from
ending = oo.seqs.xml[[1]][[feature.pos]][[i]][[interval.pos]][[1]]$GBInterval_to

qualifier.pos = match("GBFeature_quals", names(oo.seqs.xml[[1]][[feature.pos]][[i]]))

if(names(oo.seqs.xml[[1]][[feature.pos]][[i]][[qualifier.pos]])[[1]] == "GBQualifier"){
gene_name = oo.seqs.xml[[1]][[feature.pos]][[i]][[qualifier.pos]]$GBQualifier$GBQualifier_value
}
if(names(oo.seqs.xml[[1]][[feature.pos]][[i]][[qualifier.pos]])[[1]] == "$GBFeature_quals"){
gene_name = oo.seqs.xml[[1]][[feature.pos]][[i]][[qualifier.pos]]$GBFeature_quals$GBQualifier$GBQualifier_value
}

c(Type, starting, ending, gene_name)
}


})))

if(dim(Annot.seq.DF)[1] == 0){
Annot.seq.DF = as.data.frame(matrix(ncol = 4,nrow = 1))
}

# In the rare case their is no sequence and no information associated tp the accession Number
} else {

Meta.data.sequences[x,] = c(Accession.NB[x], rep(NA, 9), Accession.NB[x], NA, NA, NA)
Annot.seq.DF = as.data.frame(matrix(ncol = 4, nrow = 1))

}# end else

}# end else

names(Annot.seq.DF) = c("Type", "starting", "ending", "gene_name")
Annotation.Sequence.list[[x]] = data.frame(Access.Nb = rep(Accession.NB[x], dim(Annot.seq.DF)[1]), Annot.seq.DF)

                                            
} # end for x

colnames(Meta.data.sequences) = c("GBSeq_locus", "GBSeq_length", "GBSeq_strandedness", "GBSeq_moltype", "GBSeq_topology", "GBSeq_division", "GBSeq_update-date", "GBSeq_create-date", "GBSeq_definition", "GBSeq_primary-accession", "GBSeq_accession-version", "GBSeq_source", "GBSeq_organism", "GBSeq_taxonomy")

Annotation.Sequence.DF = do.call(rbind, Annotation.Sequence.list)

# Remove improper gene annotation for the gene.

if(length(grep("trn", Annotation.Sequence.DF$gene_name)) > 0){
Annotation.Sequence.DF = Annotation.Sequence.DF[-grep("trn", Annotation.Sequence.DF$gene_name),]
}

if(length(grep("tRNA", Annotation.Sequence.DF$gene_name)) > 0){
Annotation.Sequence.DF = Annotation.Sequence.DF[-grep("tRNA", Annotation.Sequence.DF$gene_name),]
}

if(length(grep("rrn", Annotation.Sequence.DF$gene_name)) > 0){
Annotation.Sequence.DF = Annotation.Sequence.DF[-grep("rrn", Annotation.Sequence.DF$gene_name),]
}

Mitogenome.without.gene = NA
## Report the presence of mitochondrial genome badly annotated without any gene.
if(length(which(is.na(Annotation.Sequence.DF$Type))) > 0){

Seq.withoutgene = Annotation.Sequence.DF[which(is.na(Annotation.Sequence.DF$Type)),1]
Meta.data.sequences[Meta.data.sequences[,11] %in% Seq.withoutgene,13]

Mitogenome.without.gene = data.frame(Seq.withoutgene = Seq.withoutgene, Taxa.name = Meta.data.sequences[Meta.data.sequences[,11] %in% Seq.withoutgene, 13])

} 

# Report the presence of mitochondrial genome with missing gene.
Mitogenome.with.missing.gene = NA

if(length(which(table(Annotation.Sequence.DF$Access.Nb) < 13)) > 0){

Seq.withmissinggene = labels(which(table(Annotation.Sequence.DF$Access.Nb) < 13))

Mitogenome.with.missing.gene = data.frame(Seq.with.missing.gene = Seq.withmissinggene, Taxa.name = Meta.data.sequences[Meta.data.sequences[,11] %in% Seq.withmissinggene,13])

}

# Clean the annotation
Annotation.Sequence.DF$gene_name.corrected = Annotation.Sequence.DF$gene_name

Annotation.Sequence.DF$gene_name.corrected = gsub("ox", "o", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("iii", "3", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("ii", "2", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("i", "1", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("nadh", "nd", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("nad", "nd", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub(" ", "", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("atpase", "atp", Annotation.Sequence.DF$gene_name.corrected)

Annotation.Sequence.DF$gene_name.corrected = gsub("cob", "cytb", Annotation.Sequence.DF$gene_name.corrected)


out.put = list(Meta.data.sequences = as.data.frame(Meta.data.sequences), Annotation.Sequence.DF = Annotation.Sequence.DF, Mitogenome.without.gene = Mitogenome.without.gene, Mitogenome.with.missing.gene = Mitogenome.with.missing.gene)

saveRDS(out.put, file = paste(output.Path, "/",  output.Name, ".MetaData.Annotations.rds", sep = ""))
 
return(out.put)
}

# For Debug
# rm(out.put, Mitogenome.with.missing.gene, Mitogenome.without.gene, Seq.withoutgene, Annotation.Sequence.DF, Meta.data.sequences, Annot.seq.DF, feature.pos, interval.pos, starting, ending, qualifier.pos, gene_name, Type, Align.name, oo.seqs.xml, oo.seqs, oo, Annotation.Sequence.list, output.Name, output.Path, Accession.NB)


