# Prepare GCTA matrices for reading into R
# Author: Susan Johnston

#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ set up working environment

source("r/makeGRM.R")

AnalysisSuffix <- "g"

memory.limit(size = 800000000) #ESSENTIAL OTHERWISE YOU RUN INTO MEMORY ISSUES


#~~ set working directory

setwd("gcta")

#~~ read in famped

famped <- read.table("../results/2_FamilyPedigree_afterQC_g.txt", header = T)

#~~ read the GRM file from GCTA

system.time({
  
  grm.auto <- read.table("150129_autoGRM_adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
  ids.auto <- read.table("150129_autoGRM_adj.grm.id")  # CONTAINS ID LIST
  
  ids.auto$Order <- 1:nrow(ids.auto)
  
  #~~ remove IDs not used in the analysis
  
  ids.auto <- ids.auto[which(ids.auto$V2 %in% famped[,1]),]
  grm.auto <- grm.auto[which(grm.auto$V1 %in% ids.auto$Order & grm.auto$V2 %in% ids.auto$Order),-3]
  
  save(grm.auto, ids.auto, file = "Autosomes_GRM.Rdata")
  
  rm(ids.auto, grm.auto)
  
})

system.time({
  
  grm.auto <- read.table("150129_autowparGRM_adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
  ids.auto <- read.table("150129_autoGRM_adj.grm.id")  # CONTAINS ID LIST
  
  ids.auto$Order <- 1:nrow(ids.auto)
  
  #~~ remove IDs not used in the analysis
  
  ids.auto <- ids.auto[which(ids.auto$V2 %in% famped[,1]),]
  grm.auto <- grm.auto[which(grm.auto$V1 %in% ids.auto$Order & grm.auto$V2 %in% ids.auto$Order),-3]
  
  save(grm.auto, ids.auto, file = "Autosomes_wPAR_GRM.Rdata")
  
  rm(ids.auto, grm.auto)
  
})


system.time({
  for(i in 1:26){
    grm.chr <- read.table(paste0("150129_chr", i, "_GRM_adj.grm.gz"))
    ids.chr <- read.table(paste0("150129_chr", i, "_GRM_adj.grm.id"))
    
    ids.chr$Order <- 1:nrow(ids.chr)
    
    ids.chr <- ids.chr[which(ids.chr$V2 %in% famped[,1]),]
    grm.chr <- grm.chr[which(grm.chr$V1 %in% ids.chr$Order & grm.chr$V2 %in% ids.chr$Order),-3]
    
    grm.wochr <- read.table(paste0("150129_WOchr", i, "_GRM_adj.grm.gz"))
    ids.wochr <- read.table(paste0("150129_WOchr", i, "_GRM_adj.grm.id"))
    
    ids.wochr$Order <- 1:nrow(ids.wochr)
    
    ids.wochr <- ids.wochr[which(ids.wochr$V2 %in% famped[,1]),]
    grm.wochr <- grm.wochr[which(grm.wochr$V1 %in% ids.wochr$Order & grm.wochr$V2 %in% ids.wochr$Order),-3]
    
    save(grm.chr, ids.chr, grm.wochr, ids.wochr, file = paste0("Chr_", i, "_GRMs.Rdata"))
    
    rm(grm.chr, ids.chr, grm.wochr, ids.wochr)
  }
})

#~~ sex chromosomes

# PAR

grm.par <- read.table("150129_chr27par_GRM_adj.grm.gz")
ids.par <- read.table("150129_chr27par_GRM_adj.grm.id")
ids.par$Order <- 1:nrow(ids.par)
ids.par <- ids.par[which(ids.par$V2 %in% famped[,1]),]
grm.par <- grm.par[which(grm.par$V1 %in% ids.par$Order & grm.par$V2 %in% grm.par$Order),-3]

save(grm.par, ids.par, file = paste0("Chr_27par_GRMs.Rdata"))

rm(grm.par, ids.par)

# Non-PAR

grm.nonpar <- read.table("150129_chr27nonpar_GRM.grm.gz")
ids.nonpar <- read.table("150129_chr27nonpar_GRM.grm.id")
ids.nonpar$Order <- 1:nrow(ids.nonpar)
ids.nonpar <- ids.nonpar[which(ids.nonpar$V2 %in% famped[,1]),]
grm.nonpar <- grm.nonpar[which(grm.nonpar$V1 %in% ids.nonpar$Order & grm.nonpar$V2 %in% grm.nonpar$Order),-3]

save(grm.nonpar, ids.nonpar, file = paste0("Chr_27nonpar_GRMs.Rdata"))

rm(grm.nonpar, ids.nonpar)



