

library(simperSNP)
library(reshape)
library(plyr)
library(GenABEL)

load("../DataForSimulation.Rdata")
# cmmaphold <- cmmap[sort(sample(1:39000, size = 400)),]
cmmaphold <- cmmap
rm(cmmap)

#~~ Load information file


#~~ Run Crimap



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run crimap                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ Create Crimap Files for Autosomes - do chromosome by chromosome to reduce space taken up

if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

rectab  <- NULL          # Raw data on recombination events
genabel.list <- list()

system.time({
  for(i in 1:26){
    
    message(paste("Running chromosome", i))
    
    cmmap <- subset(cmmaphold, Chr == i)
    
    system.time(chrNo <- simulateGenos(ped,  
                                       pedigree.type  = "simple",
                                       cM.male        = cmmap$cMPosition.Male,
                                       cM.female      = cmmap$cMPosition.Female,
                                       founder.mafs   = cmmap$Q.2,
                                       chromosome.ids = cmmap$Chr,
                                       snp.names      = cmmap$SNP.Name,
                                       map.distance   = cmmap$GenomicPosition,
                                       xover.min.cM   = 10))
    
    createCrimapFiles(chrNo$genabel.object,
                      chromosomeid = paste(i, AnalysisSuffix, sep = ""),
                      familyPedigree = famped,
                      #menderrtab = newerrtab,
                      chr = i)
    system(paste0("../crimap ", i, AnalysisSuffix," prepare < crimapinput1 > chr", i, AnalysisSuffix, ".pre"))
    system(paste0("../crimap ", i, AnalysisSuffix," chrompic > chr", i, AnalysisSuffix, ".cmp"))
    system("rm *.cg")
    
    chrompicfile <- paste("chr", i, AnalysisSuffix, ".cmp", sep = "")
    recombtab <- RecombRateFromChrompic(chrompicfile)
    recombtab$Chr <- i
    
    recombtab$Offspring.ID <- unlist(lapply(recombtab$Family, function(x) strsplit(x, split = "_")[[1]][3]))
    recombtab$Parent.Type <- unlist(lapply(recombtab$Family, function(x) strsplit(x, split = "_")[[1]][2]))
    recombtab$Parent.Type <- ifelse(recombtab$Parent.Type == "Mum", "MOTHER", "FATHER")
    
    recombtab <- droplevels(recombtab[which(recombtab$Offspring.ID == recombtab$id & recombtab$parent == recombtab$Parent.Type),])
    
    #~~ Add parent IDs
    transped <- melt(ped, id.vars = "ANIMAL")
    head(transped)
    names(transped) <- c("id", "parent", "Parent.ID")
    
    recombtab <- join(recombtab, transped)
    
    #~~ Get simulated recombination counts
    
    
    x <- data.frame(Template = unlist(lapply(chrNo$templates, function(x) paste(x, collapse = "_"))))
    x$Key <- names(unlist(lapply(chrNo$templates, function(x) paste(x, collapse = "_"))))
    x$Chr <- as.numeric(i)
    
    
    recombtab$Key <- paste(recombtab$Parent.ID, recombtab$Offspring.ID, sep = "_")
    recombtab <- join(recombtab, x)
    
    rectab <- rbind(rectab, recombtab)
    
    system(paste0("rm chr", i, AnalysisSuffix, "*"))
    
    genabel.list[as.numeric(i)] <- chrNo$genabel.object
  }
  
  message(paste("Running chromosome", 27))
  
  cmmap <- subset(cmmaphold, Chr == 27)
  
  system.time(chrNo <- simulateGenos(ped,  
                                     pedigree.type  = "simple",
                                     cM.male        = cmmap$cMPosition.Male,
                                     cM.female      = cmmap$cMPosition.Female,
                                     founder.mafs   = cmmap$Q.2,
                                     chromosome.ids = cmmap$Chr,
                                     snp.names      = cmmap$SNP.Name,
                                     map.distance   = cmmap$GenomicPosition,
                                     xover.min.cM   = 10))
  
  
  createCrimapFiles(chrNo$genabel.object,
                    chromosomeid = paste("X", AnalysisSuffix, sep = ""), 
                    familyPedigree = famped,
                    pseudoautoSNPs = pseudoautoSNPs,
                    chr = "X")
  
  system(paste0("mv chrX", AnalysisSuffix, ".gen chr27", AnalysisSuffix, ".gen"))
  
  system(paste0("../crimap ", 27, AnalysisSuffix," prepare < crimapinput1 > chr", 27, AnalysisSuffix, ".pre"))
  system(paste0("../crimap ", 27, AnalysisSuffix," chrompic > chr", 27, AnalysisSuffix, ".cmp"))
  system("rm *.cg")
  
  chrompicfile <- paste("chr", 27, AnalysisSuffix, ".cmp", sep = "")
  recombtab <- RecombRateFromChrompic(chrompicfile)
  recombtab$Chr <- 27
  
  recombtab$Offspring.ID <- unlist(lapply(recombtab$Family, function(x) strsplit(x, split = "_")[[1]][3]))
      recombtab$Parent.Type <- unlist(lapply(recombtab$Family, function(x) strsplit(x, split = "_")[[1]][2]))
      recombtab$Parent.Type <- ifelse(recombtab$Parent.Type == "Mum", "MOTHER", "FATHER")
      
      recombtab <- droplevels(recombtab[which(recombtab$Offspring.ID == recombtab$id & recombtab$parent == recombtab$Parent.Type),])
      
      #~~ Add parent IDs
      transped <- melt(ped, id.vars = "ANIMAL")
      head(transped)
      names(transped) <- c("id", "parent", "Parent.ID")
      
      recombtab <- join(recombtab, transped)
      
      #~~ Get simulated recombination counts
      
      
      x <- data.frame(Template = unlist(lapply(chrNo$templates, function(x) paste(x, collapse = "_"))))
      x$Key <- names(unlist(lapply(chrNo$templates, function(x) paste(x, collapse = "_"))))
      x$Chr <- as.numeric(i)
      
      
      recombtab$Key <- paste(recombtab$Parent.ID, recombtab$Offspring.ID, sep = "_")
      recombtab <- join(recombtab, x)
 
  
  
  rectab <- rbind(rectab, recombtab)
  system(paste0("rm chr", 27, AnalysisSuffix, "*"))
  
    genabel.list[as.numeric(i)] <- chrNo$genabel.object
  
  
})

#~~ Clean up data


save(rectab, genabel.list, file = paste0(AnalysisSuffix, ".Rdata"))






