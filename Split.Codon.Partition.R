#Function developed by David Eme
#' @title Add the three codons position to the partition file in nexus format 

#' @description This function use the 'Partition_Concat.nex' file provided by the 
#' Align.Concat, and provide the additional three codon positions for the coding gene in another
#' Partition_file in nexus format. This file is ready to be used by IQtree or RAxML

#' @return This function returns a Partition file in enxus format including the 
#' additional three codon positions for the coding gene
#' 
#' @param input.parti the path including the name of the Partition_file.nex exported 
#' by the Align.Concat function in nexus format.

#' @param which.codon a vector of interger providing the position of the gene in 
#' the Supermatrix that require to include the three codons positions.
#'  To know the position of the coding gene you can either look at the print stdout 
#'  from the Align.Concat function, or upload the "Convtab.txt" table provided by
#'  the Align.Concat function. 
#'  
#' @export Split.Codon.Partition

Split.Codon.Partition =  function(input.parti = NULL, which.codon = NULL, output.parti = NULL) {
  
  file = readLines(input.parti)
  Parti.ori = file[-c(1,2, length(file))]
  
  Parti.Posi.ori = unlist(lapply(strsplit(Parti.ori, " = "), function(x) x[2]))
  Parti.Posi.ori = data.frame(Seq.Order = seq(1, length(Parti.Posi.ori)), Position = Parti.Posi.ori)
  
  # insert the tree codon position
  position = match(Parti.Posi.ori[,1], which.codon)
  
  New.DF = do.call(rbind, lapply(1:dim(Parti.Posi.ori)[1], function(x){
    if(is.na(position[x])){
      res = Parti.Posi.ori[x,]
    } else {
      Start.Pos = as.numeric(strsplit(Parti.Posi.ori[which.codon[x],2], "-")[[1]][1])
      End.Pos = strsplit(Parti.Posi.ori[which.codon[x],2], "-")[[1]][2]
      
      res = cbind(Seq.Order = seq(1,3), Position = c(paste(Start.Pos, "-", gsub(";", "\\3", End.Pos, fixed= T), ";", sep = ""), paste((Start.Pos+1), "-", gsub(";", "\\3", End.Pos, fixed= T), ";", sep = ""), paste((Start.Pos+2), "-", gsub(";", "\\3", End.Pos, fixed= T), ";", sep = "")))
    }
    res
  }))
  
  New.DF = data.frame(charset_part = paste("charset part", seq(1, dim(New.DF)[1]), sep=""), Position = New.DF[,2])
  
  #Rebuild the nexus file
  cat(file = output.parti, "#nexus", "\n",
      "begin sets;", "\n", sep="")
  
  lapply(1:dim(New.DF)[1], function(x){
    cat(file = output.parti, paste(as.matrix(New.DF[x,]), collapse = " = "), "\n", sep="", append = T)
  })
  
  cat(file = output.parti, "end;", "\n", sep="", append = T)
}

