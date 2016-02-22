# Checks linkage maps against physical maps.
# Author: Susan Johnston


source("0_Chrompic_Functions.R")

library(plyr)

AnalysisSuffix <- "g"

mapdata <- NULL

for(i in 0:27) mapdata <- rbind(mapdata,cbind(MapFromChrompic(paste0("chr", i, AnalysisSuffix, ".cmp")), Chr = i))

plinkmap <- read.table("data/20130410merged1_66nodups.oar3.1.map")
plinkmap <- subset(plinkmap, select = -V3)

head(plinkmap)
names(plinkmap) <- c("Chr", "SNP.Name", "GenomicPosition")

mapdata <- join(mapdata, plinkmap)

str(mapdata)

for(i in c("Position", "r", "cMdiff")) mapdata[,i] <- as.numeric(mapdata[,i])

library(ggplot2)

ggplot(mapdata, aes(GenomicPosition, Position)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~Chr, scales = "free")

head(mapdata)
tail(mapdata)

write.table(mapdata, "results/2_Maps_from_Chrompic_g.txt", row.names = F, sep = "\t", quote = F)
