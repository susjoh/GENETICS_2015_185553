load("2.1_Analysis_data_d.Rdata")

source("0_createCrimapFilesUnix.R")

for(i in 1:26){
  
  chrNo  <- genabelSubset(sheepabel, i)
  chrNo  <- genabelSubset(sheepabel, i, snplist=snp.names(chrNo))
  createCrimapFiles(chrNo, chromosomeid = paste(i, AnalysisSuffix, sep = ""), familyPedigree=famped, directory="crimap")
  
  
  if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)
  
  if(!file.exists(paste("chr", i, AnalysisSuffix,".par", sep = ""))){
    system(paste("./crimap ", i, AnalysisSuffix, " prepare < crimapinput1", sep = ""))
  }
  
  
  #~~ Run Chrompic
  
  system(paste("./crimap ", i, AnalysisSuffix, " chrompic > chr", i, AnalysisSuffix,".cmp", sep = ""))
  
}