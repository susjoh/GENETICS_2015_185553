# Runs crimap on each chromosome. Deal with Mendelian errors.
# Author: Susan Johnston

# Requires rename in PATH and needs crimap2504.exe in the directory "crimap" (included in repo)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)


AnalysisSuffix <- "g"

source("0_Crimap_Functions_Faster_v2.R")

#~~ Load GenABEL gwaa.data for all genotyped sheep

sheepabel <- load.gwaa.data(phenofile = "data/1_GenABEL_hornpheno20150126.txt",
                            genofile =  "data/1_GenABEL_sheepabelQCed20150126.gen")

#~~ Read in pedigree file

pedigree <- read.table("data/pedigree_20130920.txt", header = T)

#~~ Read in pseudoautosomal SNP information

pseudoautoSNPs <- readLines("results/1_Pseudoautosomal_SNPs_in_X.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Determine Working Pedigrees for Crimap    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove non-genotyped parents from the pedigree

pedigree <- subset(pedigree, ANIMAL %in% idnames(sheepabel))
pedigree$MOTHER[which(!pedigree$MOTHER %in% idnames(sheepabel))] <- 0
pedigree$FATHER[which(!pedigree$FATHER %in% idnames(sheepabel))] <- 0


#~~ remove parent if there was sex-ambiguity (father or mother has different sex than from plink)

test <- phdata(sheepabel)
test$IsMum <- ifelse(test$id %in% pedigree$MOTHER, "yes", "no")
test$IsDad <- ifelse(test$id %in% pedigree$FATHER, "yes", "no")

ambiguous_sex <- c(test[test$sex == 0 & test$IsDad == "yes","id"], test[test$sex == 1 & test$IsMum == "yes","id"])

if(length(ambiguous_sex) > 0){
  pedigree <- subset(pedigree, !ANIMAL %in% ambiguous_sex)
  if(length(which(pedigree$MOTHER %in% ambiguous_sex)) > 0) pedigree$MOTHER[which(pedigree$MOTHER %in% ambiguous_sex)] <- 0
  if(length(which(pedigree$FATHER %in% ambiguous_sex)) > 0) pedigree$FATHER[which(pedigree$FATHER %in% ambiguous_sex)] <- 0
}

rm(test, ambiguous_sex)

#~~ remove parent if only one parent

pedigree$MOTHER[which(pedigree$FATHER == 0 & pedigree$MOTHER != 0)] <- 0
pedigree$FATHER[which(pedigree$FATHER != 0 & pedigree$MOTHER == 0)] <- 0


#~~ Find all offspring with two parents

offped <- pedigree[which(pedigree[,2] != 0),]
offped$MumParents <- ifelse(offped$MOTHER %in% offped$ANIMAL, "yes", "no")
offped$DadParents <- ifelse(offped$FATHER %in% offped$ANIMAL, "yes", "no")

offped <- offped[-which(offped$MumParents == "no" & offped$DadParents == "no"),]

#~~ create a Family Pedigree

famped <- NULL

for(i in offped$ANIMAL){
  
  ped1 <- offped[which(offped$ANIMAL == i),]
  
  if(ped1$MumParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree[which(pedigree$ANIMAL %in% c(ped1[,"MOTHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
      
    ped2 <- rbind(c(ped1$FATHER, 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Mum", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
  }
  
  #~~ Make founder parents
  
  if(ped1$DadParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree[which(pedigree$ANIMAL %in% c(ped1[,"FATHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    
    ped2 <- rbind(c(ped1$MOTHER, 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Dad", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
    
  }
  
  rm(ped1)
}

#~~ Remove families where an individual appears more than once

badfams <- data.frame(table(famped$ANIMAL, famped$Family))
badfams <- subset(badfams, Freq > 1)

famped <- subset(famped, !Family %in% badfams$Var2)


#~~ write family pedigree to file for future reference

write.table(famped, 
            paste("results/2_FamilyPedigree_Raw", AnalysisSuffix, ".txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)

table(table(famped$Family))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create crimap files & Run 1st instance    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!paste0("crimap_", AnalysisSuffix) %in% dir("crimap")) system(paste0("mkdir crimap_", AnalysisSuffix))
setwd(paste0("crimap/crimap_", AnalysisSuffix))

#~~ Create Crimap Files for Autosomes. NB. This is very slow (1hr) and could be parallelised, code is old but works!

system.time({
  for(i in 0:26){
    print(paste("Formatting chromosome", i))
    chrNo  <- genabelSubset(sheepabel, i)
    createCrimapFiles(chrNo, chromosomeid = paste(i, AnalysisSuffix, sep = ""), familyPedigree=famped)
  }
})

#~~ Create Crimap Files for X chromosome

print(paste("Formatting chromosome X"))
chrNo  <- genabelSubset(sheepabel, "X")
createCrimapFiles(chrNo, chromosomeid = paste("X", AnalysisSuffix, sep = ""), familyPedigree=famped, pseudoautoSNPs = pseudoautoSNPs)
system("cmd", input = paste0("rename chrX", AnalysisSuffix, ".gen chr27", AnalysisSuffix, ".gen"))

#~~ RUN CRIMAP INITIALLY TO FIND MENDELIAN INCONSISTENCIES (prepare function)

if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

for(i in c(0:27)){
  print(paste("Running chromosome", i))
  system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," prepare < crimapinput1 > chr", i, AnalysisSuffix, ".pre"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ 3. Process Mendelian inconsistencies       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

menderr <- NULL

for(i in c(0:27)){
   
  print(paste("Checking chromosome", i, "for Mendelian inconsistencies"))
  
  
  x <- readLines(paste0("chr", i, AnalysisSuffix, ".pre"))
  
  x.0 <- x[grep("NONINHERITANCE", x)]
  x.1 <- x[grep("NONINHERITANCE", x) + 1]
  x.2 <- x[grep("NONINHERITANCE", x) + 2]
  x.3 <- x[grep("NONINHERITANCE", x) + 3]
  
  if(length(x.0) > 0){
    xtab <- data.frame(Family     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][3]),
                       ANIMAL     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][5]),
                       Locus      = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][7]),
                       Mat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][7]),
                       Pat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][9]),
                       Mat.UniqueAlleles = x.2,
                       Pat.UniqueAlleles = x.3,
                       Chr       = i)
    
    row.names(xtab) <- 1:nrow(xtab)
    head(xtab)
    
    xtab$Family     <- gsub(",", "", xtab$Family)
    xtab$ANIMAL     <- gsub(",", "", xtab$ANIMAL)
    xtab$Mat.Allele <- gsub(",", "", xtab$Mat.Allele)
    xtab$Locus      <- gsub(":", "", xtab$Locus)
    xtab$Mat.UniqueAlleles <- gsub("Maternal alleles: ", "", xtab$Mat.UniqueAlleles)
    xtab$Pat.UniqueAlleles <- gsub("Paternal alleles: ", "", xtab$Pat.UniqueAlleles)
    
    
    menderr <- rbind(menderr, xtab)
  }
  
  rm(x.0, x.1, x.2, x.3, xtab)
}


#~~ Merge with pedigree information

menderr <- join(menderr, pedigree)

#~~ Identify individuals with many mismatches

menderr$Mat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Mat.Allele, menderr$Mat.UniqueAlleles)
menderr$Pat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Pat.Allele, menderr$Pat.UniqueAlleles)

parent.mismatch <- data.frame(table(menderr$ANIMAL, menderr$Mat.Mismatch, menderr$Pat.Mismatch))

parent.mismatch <- unique(parent.mismatch[-which(parent.mismatch$Freq == 0),])

head(parent.mismatch)
names(parent.mismatch) <- c("ANIMAL", "Mat.Mismatches", "Pat.Mismatches", "Count")

test <- rbind(cbind(parent.mismatch[which(parent.mismatch$Mat.Mismatches == "yes"), c("ANIMAL", "Count")],
                    Parent = "MOTHER"),
              cbind(parent.mismatch[which(parent.mismatch$Pat.Mismatches == "yes"), c("ANIMAL", "Count")],
                    Parent = "FATHER"))

head(test)

parent.mismatch <- test

#~~ As there are repeated IDs in there, extract the highest incidence of the count

parent.mismatch$Parent.ID.Code <- paste(parent.mismatch$ANIMAL, parent.mismatch$Parent, sep = "_")
parent.dup.tab <- data.frame(table(parent.mismatch$Parent.ID.Code))
parent.dup.tab <- subset(parent.dup.tab, Freq > 1)

for(i in 1:nrow(parent.dup.tab)){
  count.temp <- parent.mismatch$Count[which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i])]
  
  if(count.temp[1] != count.temp[2]){
    parent.mismatch <- parent.mismatch[-which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i] & 
                                                parent.mismatch$Count == min(count.temp)),]
  }
  
  if(count.temp[1] == count.temp[2]) {
    
    parent.mismatch <- parent.mismatch[-which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i] & 
                                                parent.mismatch$Count == min(count.temp))[1],]
  }
}  


exclusion.threshold <- nsnps(sheepabel) * 0.001

ggplot(parent.mismatch[which(parent.mismatch$Count < exclusion.threshold),], aes(Count)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~Parent)

remove.ids <- parent.mismatch[which(parent.mismatch$Count > exclusion.threshold),]

remove.ids$ANIMAL <- as.numeric(as.character(remove.ids$ANIMAL))
remove.ids$Parent <- as.character(remove.ids$Parent)


#~~ revise the pedigree and error table

pedigree.new <- pedigree

if(nrow(remove.ids) > 0){
  for(i in 1:nrow(remove.ids)){
    
    pedigree.new[which(pedigree.new$ANIMAL == remove.ids$ANIMAL[i]), remove.ids$Parent[i]] <- 0
    
    if(remove.ids$Parent[i] == "MOTHER") menderr <- menderr[-which(menderr$ANIMAL == remove.ids$ANIMAL[i] & menderr$Mat.Mismatch == "yes"),]
    if(remove.ids$Parent[i] == "FATHER") menderr <- menderr[-which(menderr$ANIMAL == remove.ids$ANIMAL[i] & menderr$Pat.Mismatch == "yes"),]
        
  }
}

#~~ Deal with problematic loci within parents by removing genotype in parents and offspring.

head(menderr)

newerrtab <- menderr[,c("ANIMAL", "Locus", "Chr")]

mattab <- menderr[which(menderr$Mat.Mismatch == "yes"),c("MOTHER", "Locus", "Chr")]
names(mattab) <- c("ANIMAL", "Locus", "Chr")
pattab <- menderr[which(menderr$Pat.Mismatch == "yes"),c("FATHER", "Locus", "Chr")]
names(pattab) <- c("ANIMAL", "Locus", "Chr")

newerrtab <- rbind(newerrtab, mattab, pattab)

tapply(newerrtab$Locus, newerrtab$Chr, max)

head(newerrtab)

newerrtab$SNP.Name <- NA

for(i in 1:nrow(newerrtab)) newerrtab$SNP.Name[i] <- names(which(chromosome(sheepabel) == newerrtab$Chr[i])[newerrtab$Locus[i]])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Determine Working Pedigrees for Crimap    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove non-genotyped parents from the pedigree

pedigree.new$MOTHER[which(!pedigree.new$MOTHER %in% idnames(sheepabel))] <- 0
pedigree.new$FATHER[which(!pedigree.new$FATHER %in% idnames(sheepabel))] <- 0


#~~ remove parent if only one parent
pedigree.new$MOTHER[which(pedigree.new$FATHER == 0 & pedigree.new$MOTHER != 0)] <- 0
pedigree.new$FATHER[which(pedigree.new$FATHER != 0 & pedigree.new$MOTHER == 0)] <- 0

#~~ Find all offspring with two parents

offped <- pedigree.new[which(pedigree.new[,2] != 0),]
offped$MumParents <- ifelse(offped$MOTHER %in% offped$ANIMAL, "yes", "no")
offped$DadParents <- ifelse(offped$FATHER %in% offped$ANIMAL, "yes", "no")

#~~ create a Family pedigree.new

famped <- NULL

for(i in offped$ANIMAL){
  
  ped1 <- offped[which(offped$ANIMAL == i),]
  
  if(ped1$MumParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree.new[which(pedigree.new$ANIMAL %in% c(ped1[,"MOTHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    
    ped2 <- rbind(c(ped1$FATHER, 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Mum", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
  }
  
  #~~ Make founder parents
  
  if(ped1$DadParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree.new[which(pedigree.new$ANIMAL %in% c(ped1[,"FATHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    
    ped2 <- rbind(c(ped1$MOTHER, 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Dad", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
    
  }
  
  rm(ped1)
}

#~~ write family pedigree to file for future reference

famped <- subset(famped, !Family %in% badfams$Var2)


write.table(famped, paste0("../../results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(newerrtab, paste0("../../results/2_MendelianErrors_Crimap_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Rerun crimap                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("cmd", input = "del chr*")

#~~ Create Crimap Files for Autosomes

system.time({
  for(i in 0:26){
    chrNo  <- genabelSubset(sheepabel, i)
    createCrimapFiles(chrNo,
                      chromosomeid = paste(i, AnalysisSuffix, sep = ""),
                      familyPedigree = famped,
                      menderrtab = newerrtab,
                      chr = i)
  }
})

#~~ Create Crimap Files for X chromosome

chrNo  <- genabelSubset(sheepabel, "X")
createCrimapFiles(chrNo,
                  chromosomeid = paste("X", AnalysisSuffix, sep = ""), 
                  familyPedigree = famped,
                  pseudoautoSNPs = pseudoautoSNPs,
                  menderrtab = newerrtab,
                  chr = 27)

system("cmd", input = paste0("rename chrX", AnalysisSuffix, ".gen chr27", AnalysisSuffix, ".gen"))

if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

for(i in 0:27){
  system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," prepare < crimapinput1 > chr", i, AnalysisSuffix, ".pre"))
}


for(i in 0:27){
  if(length(grep("NONINHERITANCE", readLines(paste0("chr", i, AnalysisSuffix, ".pre")))) > 0) stop()
}



for(i in 0:27){
  system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," chrompic > chr", i, AnalysisSuffix, ".cmp"))
  system("cmd", input = "del *.cg")
}


