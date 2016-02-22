# Extract the sex-averaged and sex-specific linkage maps using crimap
# Author: Susan Johnston

# Requires crimap2504.exe in crimap directory

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(plyr)

AnalysisSuffix <- "g"

setwd(paste0("crimap/crimap_", AnalysisSuffix))
source("../../0_Chrompic_Functions.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run CRIMAP to generate .map files         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for(i in c(1:27)){
 
  print(paste("Running chromosome", i))
   
  N <- read.table(paste0("chr", i, AnalysisSuffix, ".gen"), skip = 1, nrows = 1)[[1]]
  
  write.table(1, paste0("chr", i, AnalysisSuffix, ".ord"), quote = F, row.names = F, col.names = F)
  write.table("", paste0("chr", i, AnalysisSuffix, ".ord"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste(c(1, N), collapse = " "), paste0("chr", i, AnalysisSuffix, ".ord"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste(c(0:(N - 1)), collapse = " "), paste0("chr", i, AnalysisSuffix, ".ord"), quote = F, row.names = F, col.names = F, append = T)
  
  print(system.time(system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," build > chr", i, AnalysisSuffix, ".map"))))
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Extract map information from .map and .cmp files         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fullmap <- NULL

for(i in 1:27){
  
  #for(i in c(24)){
  
  tempmap <- MapFromMap(paste0("chr", i, AnalysisSuffix, ".map"))
  tempmap$Chr <- i
  
  recombmap <- MapFromChrompic(paste0("chr", i, AnalysisSuffix, ".cmp"))
  recombmap$Chr <- i
  
  recombmap <- join(recombmap, tempmap[,-which(names(tempmap) %in% c("Order", "Chr"))])
  
  fullmap <- rbind(fullmap, recombmap)
  rm(tempmap, recombmap)
}

names(fullmap)[which(names(fullmap) == "Position")] <- "cMPosition"


fullmap[,which(names(fullmap) != "SNP.Name")] <- lapply(fullmap[,which(names(fullmap) != "SNP.Name")] , as.numeric)



#~~ Deal with the X-sex-averaged map

head(fullmap)

fullmap$new.cM <- apply(fullmap[,c("cMPosition.Male", "cMPosition.Female")], 1, mean)

#~~ write to file

write.table(fullmap, paste0("../../results/2_Linkage_Map_Positions_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")

library(ggplot2)

snpmap <- read.table("../../data/20130410merged1_66nodups.oar3.1.map")
names(snpmap) <- c("Chr", "SNP.Name", "X", "Position")
snpmap <- subset(snpmap, select = -X)

library(plyr)
library(reshape)

plotmap <- fullmap[,c("SNP.Name", "cMPosition.Female", "cMPosition.Male", "new.cM")]

names(plotmap)[2:4] <- c("Female", "Male", "Both")


plotmap <- melt(plotmap, id.vars = c("SNP.Name"))
names(plotmap)[2:3] <- c("SEX", "cM.Position")

plotmap <- join(plotmap, snpmap)

head(plotmap)
plotmap$cM.Position <- as.numeric(as.character(plotmap$cM.Position))

ggplot(plotmap, aes(Position, cM.Position, col = SEX)) +
  geom_point(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Chr, scales = "free")

Xmap <- subset(map, Chr == 27)

ggplot(Xmap, aes(Position, cM.Position, col = SEX)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Chr, scales = "free")
