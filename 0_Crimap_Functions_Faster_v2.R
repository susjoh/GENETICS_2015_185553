# Functions for creating crimap files from GenABEL objects
# Author: Susan Johnston

# Prerequisites (need to be introduced for data checking)
# pedigree must be in format ANIMAL FATHER MOTHER with those headers
# unix tools required on system: fgrep, sed


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. BASIC FUNCTIONS                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

corner <- function(data.frame) data.frame[1:5, 1:5]

save.test <- function(object) write.table(object, "test.txt", quote = F, sep = "\t", row.names = F)

#~~ GGPLOT2 default options

library(ggplot2)

defopts <- theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 strip.background = element_blank())

#~~ recoder func

recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ genabelSubset: Make a chromosome object with pedigree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FUNCTION: Make a subset of a gwaa.data object by chromosome and ID.

# INPUT: 
# gwaa.data: GenABEL gwaa.data object
# chr      : Chromosome Number(s) as vector
# pedigree : data.frame with pedigree with column names c("ANIMAL", "FATHER", "MOTHER")
# sampleno : (default = NULL) number of SNP loci to sample from the dataset
# ordsample: (default = TRUE) whether to order sampled SNP loci by their position on the chromosome

# OUTPUT:
# chrobj: GenABEL gwaa.data object


genabelSubset <- function(gwaa.data, chr = NULL, pedigree = NULL, sampleno = NULL, snplist = NULL, ordsample = TRUE){
  
  require(GenABEL)
  
  if(!is.null(chr)){
    
    chrobj <- gwaa.data[,which(chromosome(gwaa.data) %in% chr)]    #~~ Subset the data by Chromosome
    
  }
  
  if(!is.null(pedigree)){
    
    chrobj <- chrobj[which(idnames(chrobj) %in% pedigree[,1]),]  #~~ subset the data further by IDs in the pedigree object
    
  }
  
  #~~ extract snps if list is specified
  
  if(!is.null(snplist)){
    
    chrobj <- chrobj[,snplist]
    
  }
  
  #~~ conduct SNP sampling if specified
  
  if(!is.null(sampleno)){
    
    subsnps <- sample(1:nsnps(chrobj), size = sampleno, replace=F)
    
    if(ordsample == TRUE) subsnps <- sort(subsnps)
    
    chrobj <- chrobj[,subsnps]
    
  }
  
  
  
  #~~ return the new gwaa.data object
  
  chrobj
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# genabel2PreCrimap: Create a data frame in pre-Crimap format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


genabel2PreCrimap <- function(gwaa.data, pedigree){
  
  pedigree <- pedigree[,1:3]
  # FUNCTION: Create a data frame in pre-Crimap format
  
  # INPUT: 
  # gwaa.data: GenABEL gwaa.data object
  # pedigree : data.frame with pedigree with column names c("ANIMAL", "FATHER", "MOTHER")
  
  
  # OUTPUT:
  # test: data.frame in pre-Crimap format
  
  require(GenABEL)
  require(data.table)
  
  #~~ select only individuals which are within the pedigree provided
  
  ids <- as.character(unique(c(pedigree[,1], pedigree[,2], pedigree[,3])))
  ids <- ids[which(ids %in% idnames(gwaa.data))]
  
  gwaa.data <- gwaa.data[ids,]
  
  #~~ Recode ACGT as 1234 with space delimited between alleles
  
  test <- as.character.snp.data(gtdata(gwaa.data))
  test <- as.matrix(test)
  
  print("Recoding alleles to numeric values...")
  
  test[which(is.na(test))] <- "0 0"
  test[which(test == "A/A")] <- "1 1"
  test[which(test == "A/C")] <- "1 2"
  test[which(test == "A/G")] <- "1 3"
  test[which(test == "A/T")] <- "1 4"
  test[which(test == "C/A")] <- "2 1"
  test[which(test == "C/C")] <- "2 2"
  test[which(test == "C/G")] <- "2 3"
  test[which(test == "C/T")] <- "2 4"
  test[which(test == "G/A")] <- "3 1"
  test[which(test == "G/C")] <- "3 2"
  test[which(test == "G/G")] <- "3 3"
  test[which(test == "G/T")] <- "3 4"
  test[which(test == "T/A")] <- "4 1"
  test[which(test == "T/C")] <- "4 2"
  test[which(test == "T/G")] <- "4 3"
  test[which(test == "T/T")] <- "4 4"
  test[which(test == "1/1")] <- "1 1"
  test[which(test == "1/2")] <- "1 2"
  test[which(test == "2/1")] <- "2 1"
  test[which(test == "2/2")] <- "2 2"
  
  
  
  print("...done.")
  
  
  #~~ Create a data frame with ID, FATHER, MOTHER, SEX, and Genotypes.
  
  test <- data.frame(test, stringsAsFactors=F)
  
  test$ANIMAL <- idnames(gwaa.data)
  test$SEX <- phdata(gwaa.data)$sex
  
  
  pedigree$ANIMAL <- as.character(pedigree$ANIMAL)
  test2 <- data.table(test, key = "ANIMAL")
  ped2 <- data.table(pedigree, key = "ANIMAL")
  
  
  print("Merging pedigree and genotype information...")
  system.time(test2 <- merge(test2, ped2, all.y = T))
  print("...done.")
  
  test2 <- data.frame(test2)
  test2 <- unique(test2)

  
  system.time(test3 <- cbind(test2[,c("ANIMAL", "MOTHER", "FATHER", "SEX")], test2[,which(!names(test2) %in% c("ANIMAL", "MOTHER", "FATHER", "SEX"))]))
  
  #~~ put in missing sexes
  
#   test3$SEX[which(test3$ANIMAL %in% pedigree$MOTHER)] <- 0
#   test3$SEX[which(test3$ANIMAL %in% pedigree$FATHER)] <- 1
#   test3$SEX[which(is.na(test3$SEX))] <- 3
#   
  #~~ Deal with NA's
  
  corner(test3)
  test3$MOTHER[which(is.na(test3$MOTHER))] <- 0
  test3$FATHER[which(is.na(test3$FATHER))] <- 0
    
  for(i in 1:ncol(test3)) test3[which(is.na(test3[,i])),i] <- "0 0"
  
  #~~ return object
  
  test3
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create Crimap Files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# gwaa.data <- chrNo
# chromosomeid <- paste(i, AnalysisSuffix, sep = "")
# familyPedigree <- famped
# # menderrtab <- newerrtab
# chr = 24

createCrimapFiles <- function(gwaa.data, chromosomeid, familyPedigree, pseudoautoSNPs, menderrtab = NULL, chr = NULL, outdir = NULL) {
  
  nfamilies <- length(unique(familyPedigree$Family))  
  nloci <- nsnps(gwaa.data)
  locus.names <- snp.names(gwaa.data)
  
  if(!is.null(outdir))  outfile <- paste0(outdir, "/chr", chromosomeid, ".gen") else outfile <- paste0("chr", chromosomeid, ".gen")
  
  write.table(nfamilies, outfile, row.names = F, quote = F, col.names = F)
  write.table(nloci, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(locus.names, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
  
  
  #~~ create a master file of genotypes
  system.time(genotab <- genabel2PreCrimap(gwaa.data, pedigree=familyPedigree[,1:3]))

  genotab <- unique(genotab)
  
  genotab <- merge(familyPedigree, genotab, all.x = T)
  
  genotab <- genotab[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", names(genotab)[6:ncol(genotab)])]


  if(length(grep("^X", chromosomeid)) > 0){
    xmarkers <- which(!names(genotab) %in% pseudoautoSNPs)
    xmarkers <- xmarkers[6:length(xmarkers)]
    
    heterozygotes <- c("1 2", "1 3", "1 4", "2 1", "2 3", "2 4", "3 1", "3 2", "3 4", "4 1", "4 2", "4 3")

    for(j in xmarkers){
      
      #~~ remove heterozygotes
      
      genotab[which(genotab$SEX == 1 & genotab[,j] %in% heterozygotes),j] <- "0 0"
      
      genotab[which(genotab$SEX == 1 & genotab[,j] == "1 1"),j] <- "1 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "2 2"),j] <- "2 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "3 3"),j] <- "3 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "4 4"),j] <- "4 0"
      
    }
    
    rm(xmarkers, heterozygotes)
    
  }
  
    
  #~~ deal with mendelian errors
  
  if(!is.null(menderrtab)){
    menderrtab <- subset(menderrtab, Chr == chr)
    if(nrow(menderrtab) > 0){
      for(i in 1:nrow(menderrtab)){
        genotab[which(genotab$ANIMAL == menderrtab$ANIMAL[i]), 6 + as.numeric(menderrtab$Locus[i])] <- "0 0"
      }
    }
  }
          
  genotab <- as.matrix(genotab)
          
  
  counter <- 0
  
  for(i in unique(familyPedigree$Family)){
    
    counter <- counter + 1
    if(counter %in% seq(1, length(unique(familyPedigree$Family)), 200)) {
      print(paste("Analysing Family", counter, "of", length(unique(familyPedigree$Family))))
    }

    famtab <- genotab[which(genotab[,"Family"] == i),-1]
    famsize <- nrow(famtab)
    
    write.table(i, outfile, row.names = F, quote = F, col.names = F, append=T)
    write.table(famsize, outfile, row.names = F, quote = F, col.names = F, append=T)
    write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
    write.table(famtab, outfile, row.names = F, quote = F, col.names = F, append=T)
    write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
    
  }

  
}


