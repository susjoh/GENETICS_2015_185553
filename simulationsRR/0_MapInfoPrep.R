
library(GenABEL)
library(plyr)

sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")
ped <- read.table("data/pedigree_20130920_13.txt", header = T)
cMmap <- read.table("results/2_Merged_map_g.txt", header = T)
cMmap <- subset(cMmap, select = c(SNP.Name, Chr, cMPosition.Male, cMPosition.Female, GenomicPosition))

maf.info <- summary.snp.data(gtdata(sheepabel))
head(maf.info)
maf.info$SNP.Name <- row.names(maf.info)

maf.info <- subset(maf.info, select = c(SNP.Name, Q.2))

cMmap <- join(cMmap, maf.info)
head(cMmap)


write.table(cMmap, "simulationsRR/1_MapDataForRRsimulation.txt", row.names = F, sep = "\t", quote = F)


source("0_Crimap_Functions_Faster_v2.R")
source("0_Chrompic_Functions.R")
source("0_RecCount_Informative_Length_Funcs.R")
cmmap <- read.table("simulationsRR/1_MapDataForRRsimulation.txt", header = T)
ped   <- read.table("data/pedigree_20130920.txt", header = T)

pseudoautoSNPs <- readLines("results/1_Pseudoautosomal_SNPs_in_X.txt")
famped <- read.table("results/2_FamilyPedigree_afterQC_g.txt", header = T)

rm(sheepabel, cMmap, maf.info, defopts, save.test)

save(list = ls(all.names = F), file = "simulationsRR/DataForSimulation.Rdata")

ls(all.names = F)
