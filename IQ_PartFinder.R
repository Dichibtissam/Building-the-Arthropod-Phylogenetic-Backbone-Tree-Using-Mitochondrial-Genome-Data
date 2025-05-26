#' @title Estimate the best number of partitions and substitution model for each partition: wrapper for PartitionFinder implemented in IQtree2

#' @description This function is a wrapper running PartitionFinder (Lanfear et al. 2017) implemented in IQtree2 to estimate the best number of partitions and their best associated substitution models.

#'@details IQtree2 http://www.iqtree.org/ must be install and in the path and could be run using "iqtree2"
#'
#' Information about the different options to run partitionFinder in IQtree2.
#' -s "input.al" alignment
#' -p "input.part" partition file following the nexus format or the RAxML format see http://www.iqtree.org/doc/Complex-Models#partition-models
#' -m option for modelFinder, options 1)"TESTONLY" = Standard model selection (like jModelTest, ProtTest), 2) "TEST" = Standard model selection followed by tree inference, 3) "MF" =  Standard model selection with FreeRate heterogeneity (Heterotachy), 4) "MFP" = Extended model selection followed by tree inference, 5) "MFP+MERGE" = Extended model selection followed by merging the different partition until the model fit does not increase any further (best partition scheme) and then perform the tree inference, 6) "MF+MERGE" = Standard model selection with FreeRate heterogeneity follow by the selection of the best partition scheme, 7) "TESTMERGE" = Standard model selection (like jModelTest, ProtTest), follow by the selection of the best partition scheme, and the tree inference, 8) "TESTMERGEONLY" = Standard model selection (like jModelTest, ProtTest), follow by the selection of the best partition scheme (no tree inference, so equivalent of standar PartitionFinder2 analysis), 9) a string with any substitution model HKY (default), JC, F81, K2P, .... GTR, model can be completed with state frequency and rate heterogeneity (+I or +G or + R)... see iqtree2 --help for the list of all option or see page http://www.iqtree.org/doc/Substitution-Models
#'
#' --mset Restrict search to models supported by other programs (raxml, phyml, mrbayes, beast1 or beast2)
#' -T "Threads" allow to run IQtree in parallel, AUTO = automatic selection of the number of thread, or a number with a specified number of threads
#' -ntmax the maximum number of thread to use with the AUTO option.
#' -merge greedy|rcluster|rclusterf Set merging algorithm (default: rclusterf)
#' --merit AIC|AICc|BIC  Akaike|Bayesian information criterion (default: BIC)


#' @param path.input the name of the path of the input folder (but also the output folder).

#' @param input.al the name (without the path) to the input alignments in fasta format.

#' @param the name (without the path) of the file providing the gene partition of the alignment file. This file follows the nexus format or the RAxML format see http://www.iqtree.org/doc/Complex-Models#partition-models

#' @param m option for modelFinder, can be TESTONLY, TEST, MF, MFP, MFP+MERGE, MF+MERGE, MF+MERGE, TESTMERGE, TESTMERGEONLY, or a string  with any substitution model HKY (default), JC, F81, K2P, .... GTR, model can be completed with state frequency and rate heterogeneity (+I or +G or + R)... see iqtree2. See Details, default "TESTMERGEONLY",

#' @param mset Restrict search to models supported by other programs (raxml, phyml, mrbayes, beast1 or beast2) (default beast2).

#' @param Threads allow to run IQtree in parallel, AUTO = automatic selection of the number of thread, or a number with a specified number of threads, default "AUTO".

#' @param ntmax the maximum number of thread to use with the AUTO option, (default 4).

#' @param merge greedy|rcluster|rclusterf Set merging algorithm (default: rclusterf)

#' @param merit AIC|AICc|BIC  Akaike|Bayesian information criterion (default: BIC)

#' @param prepare.files if TRUE then it prepare alignment files in nexus format corresponding to the best partitioning Scheme and a file with the best partitionning substitution model for each partition. It also remove all the intermediate files.

#' @param clean if TRUE (default), the function auatomatically clean (removed) the intremediate files.

#' @return the function returns the alignments in nexus format for 

IQ_PartFinder = function(path.input = NULL, input.al = NULL, input.part = NULL, m = "TESTMERGEONLY", mset = "beast2", Threads = "AUTO", ntmax = 4, merge = "rclusterf", 
merit = "BIC", prepare.files = TRUE, clean = TRUE){

### prepare 
if(is.null(mset)){
comd.iqtree = paste("iqtree2 -s ",  input.al, " -p ", input.part, " -m ", m, " -T ", Threads, " -ntmax ", ntmax, " --merge ", merge, " --merit ", merit, sep="")
} else {
comd.iqtree = paste("iqtree2 -s ", input.al, " -p ", input.part, " -m ", m, " --mset ", mset, " -T ", Threads, " -ntmax ", ntmax, " --merge ", merge, " --merit ", merit, sep="")
}


cat(file = "Test.sh", "#!/bin/bash", "\n", "\n", "\n", comd.iqtree, "\n", sep = "")
system("chmod u+x Test.sh")

file.copy("Test.sh", path.input, overwrite = TRUE)
system(paste("cd ", path.input, "\n","./Test.sh", sep = ""))
file.remove("Test.sh")


if(prepare.files == TRUE){

PartScheme = readLines(paste(path.input, "/", input.part, ".best_model.nex", sep = ""))
Parti.nb = PartScheme[grep("charset", PartScheme)]

# Load the iniital alignment
#Align = read.dna(paste(path.input, "/", input.al, sep = ""), format = "fasta") 
Align = seqinr::read.fasta(paste(path.input, "/", input.al, sep = ""))
Align = ape::as.alignment(Align)


# convert the sequence into a matrix with every columns a 
dna.mat = do.call(rbind, lapply(1:length(Align$seq), function(i){
seqinr::s2c(Align$seq[i])
}))



lapply(1:length(Parti.nb), function(x){
pos.number = strsplit(strsplit(Parti.nb[x], " = , ", fixed = T)[[1]][2], "  ", fixed = T)[[1]]

Pos.2.keep = unlist(lapply(1:length(pos.number), function(i){

Codon.part = grep("\\", pos.number[i], fixed = T)
start = as.numeric(strsplit(pos.number[i],"-", fixed = T)[[1]][1])

if(length(Codon.part) > 0){
end = as.numeric(gsub("\\3", "",  strsplit(gsub(";", "", pos.number[i]),"-", fixed = T)[[1]][2], fixed = T))
seq(start, end, 3)
} else {
end = as.numeric(gsub(";", "", strsplit(pos.number[i],"-", fixed = T)[[1]][2], fixed = T))
seq(start, end)
}
}))

### Select the position 
Select.pos.DNA = dna.mat[,Pos.2.keep]


New.align = list()
New.align$nb = Align$nb
New.align$seq = unlist(lapply(1:dim(Select.pos.DNA)[1], function(i){
seqinr::c2s(Select.pos.DNA[i,])
}))
#New.align$nam = Align$nam
#New.align$com = NA
#class(New.align) = "alignment"

SuperMatDF = cbind(Seq.name = Align$nam, Sequences = New.align$seq)

ConcatName = paste(path.input, "/", input.al, ".part", x, ".nex", sep = "")

NBChar = nchar(New.align$seq[1])

    cat(file = ConcatName, "#NEXUS", "\n", "\n", "BEGIN DATA;", "\n", "\t", paste("DIMENSIONS NTAX=", New.align$nb, sep = ""), paste(" NCHAR=", NBChar, ";", sep = ""), sep = "",
        append = TRUE)
    cat(file = ConcatName, "\n", "\t", "FORMAT DATATYPE=DNA GAP=-;", "\n", "MATRIX",
        "\n", sep = "", append = TRUE)
    utils::write.table(SuperMatDF, file = ConcatName, sep = "\t", append = TRUE, col.names = FALSE,
        row.names = FALSE, quote = FALSE)
    cat(file = ConcatName, "\t", ";", "\n", "END;", sep = "", append = TRUE)

})
} # end prepare.files

# option to clean the output folder and remove the intermediate files
if(clean == TRUE){
unlink(paste(path.input, "/", input.part, ".ckp.gz", sep = ""))
unlink(paste(path.input, "/", input.part, ".uniqueseq.phy", sep = ""))
unlink(paste(path.input, "/", input.part, ".best_model.nex", sep = ""))
unlink(paste(path.input, "/", input.part, ".best_scheme", sep = ""))
unlink(paste(path.input, "/", input.part, ".log", sep = ""))
unlink(paste(path.input, "/", input.part, ".iqtree", sep = ""))
unlink(paste(path.input, "/", input.part, ".treefile", sep = ""))
unlink(paste(path.input, "/Test.sh", sep = ""))
unlink(paste(path.input,  "/", input.part, ".model.gz", sep = ""))
}# end clean

}

#' Help for debug
#' rm(path.input, input.al, input.part, m, mset, Threads, ntmax, merge, merit, ConcatName, NBChar, Select.pos.DNA, New.align, pos.number, dna.mat, Pos.2.keep, PartScheme, Parti.nb)




