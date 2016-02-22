
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ set working directory


#~~ load functions and libraries

library(beepr)
library(ggplot2)
library(asreml)
library(GenABEL)
library(reshape)
library(plyr)


options(stringsAsFactors = F)

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in data

# Individual summed data
recsumm  <- read.table(paste0("results/2_TotalSummIndivRR_FullCleanPostSim_", AnalysisSuffix, ".txt"), header = T, sep = "\t")

factorvec <- c("RRID", "Offspring.ID", "UniqueID2", "Family", "RRID.SEX", 
               "RRID.BYEAR", "Offspring.SEX", "Offspring.BYEAR", 
               "RRID.CAPAGE", "MateID")

recsumm[,factorvec] <- lapply(recsumm[,factorvec], as.factor)
rm(factorvec)

rectab     <- read.table(paste0("results/2_IndivRR_FullClean_", AnalysisSuffix, ".txt"), header = T, sep = "\t")

# rectab
factorvec <- c("Chr", "Offspring.ID", "parent", "Family", "RRID", "FATHER",
               "MOTHER", "RRID.SEX", "RRID.BYEAR", "Offspring.SEX",
               "Offspring.BYEAR", "RRID.CAPAGE", "RRID.Offspring.ID",
               "UniqueID", "UniqueID2", "MateID")

rectab[,factorvec] <- lapply(rectab[,factorvec], as.factor)
rm(factorvec)

names(rectab)

rectab <- rectab[,c("Family", "RRID", "RRID.SEX", "RRID.BYEAR", "RRID.CAPAGE", "RRID.Fhat3", 
                    "MateID", "Mate.Fhat3", "Offspring.ID", "Offspring.SEX", "Offspring.BYEAR", 
                    "Offspring.Fhat3", "Mean.FamR.wids", 
                    "RRID.Offspring.ID", "UniqueID", "UniqueID2", "Chr", "Chromosome.Length", 
                    "Chromosome.SNP.Count", "No.Inf.Loci", "Inf.Chr.Length", "First.Inf.Pos", 
                    "Last.Inf.Pos", "Prop.Inf.Chromosome", "Prop.Inf.SNPs", "data.v2", 
                    "RecombCount.v2", "RecombRate")]

rectab <- subset(rectab, select = c(RRID, UniqueID2, RRID.SEX, RRID.Fhat3, Chr, RecombCount.v2))


# Pedigree
pedigree <- read.table("data/pedigree_20130920.txt", header = T)

# GenABEL data
sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")

table(chromosome(sheepabel))

sheepabel <- recodeChromosome(sheepabel, list("27"="X"))
nsnps(sheepabel)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Prepare input for Analysis                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ format pedigree

names(pedigree)[1] <- "RRID"
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])
dim(pedigree)

#~~ Get a list of analysis f ids


memory.limit(size = 800000000) #ESSENTIAL OTHERWISE YOU RUN INTO MEMORY ISSUES

ainv <- asreml.Ainverse(pedigree)$ginv

# load("gcta/Autosomes_GRM.Rdata")

# ALL IDS
newrecsumm <- dplyr::select(recsumm, TotalRecombCount, RRID.SEX, RRID.Fhat3, RRID)


save(ainv, newrecsumm, pedigree, recsumm, sheepabel, rectab, file = "mixedmodelGWAS/GLMM.data.Rdata")
