# Extract information from the chrompic files
# Author: Susan Johnston

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

AnalysisSuffix <- "g"


#~~ set working directory

setwd(paste0("crimap/crimap_", AnalysisSuffix))

#~~ load functions and libraries

library(beepr)
library(ggplot2)
library(plyr)

source("../../0_Chrompic_Functions.R")
source("../../0_RecCount_Informative_Length_Funcs.R")
source("r/recoderFunc.R")

#~~ set analysis parameters


#~~ read in genomic map positions

snpmap <- read.table("../../data/20130410merged1_66nodups.oar3.1.map")
snpmap <- snpmap[,c(1, 2, 4)]
names(snpmap) <- c("Chr", "SNP.Name", "GenomicPosition")

#~~ read in linkage map positions

maptab <- read.table(paste0("../../results/2_Linkage_Map_Positions_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ read in pedigree

crimapped <- read.table(paste0("../../results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T)

#~~ read in genomic inbreeding coefficients

inbredtab <- read.table("../../gcta/20150129_autosomal_IBD.ibc", header = T)

#~~ add basic data

basedata <- read.table("../../data/20140210_SoayBaseData.txt", header = T, sep = "\t")
names(basedata) <- toupper(names(basedata))
head(basedata)

basedata$SEX <- recoderFunc(basedata$SEX, c(1, 2), c("Female", "Male"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract results from Chrompic Files       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create tables in which to store results

rectab  <- NULL          # Raw data on recombination events
doubtab <- NULL          # Information on double recombinants

#~~ Extract linkage maps from the chrompic files

system.time({
 for(i in c(1:27)){
  
    print(paste("Running Chromosome", i))
    
    chrompicfile <- paste("chr", i, AnalysisSuffix, ".cmp", sep = "")
    doublefile   <- paste("chr", i, AnalysisSuffix, ".double", sep = "")
    
    #~~ get recombination information
    
    recombtab <- RecombRateFromChrompic(chrompicfile)
    recombtab$Chr <- i
    rectab <- rbind(rectab, recombtab)
    
    
    #~~ read in information on double recombinants
    
    doubletab <- read.table(doublefile, header = T)
    doubletab$Chr <- i
    doubtab <- rbind(doubtab, doubletab)
  }
})


beep()

write.table(doubtab,
            paste("../../results/2_Double_Rec_Tab_", AnalysisSuffix, ".txt", sep = ""), 
            sep = "\t", quote = F, row.names = F)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Merge the linkage and genomic maps together #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lmap <- merge(maptab, snpmap, all.x=T)

#~~ determine the distance between each locus for genomic position and in cM

lmap <- lmap[with(lmap, order(Chr, new.cM, GenomicPosition)),]

lmap$GenomicDiff <- c(diff(lmap$GenomicPosition), 0)
lmap$cMDiff      <- c(diff(lmap$new.cM),        0)

lmap$GenomicDiff[which(lmap$GenomicDiff < 0)] <- 0
lmap$cMDiff[which(lmap$cMDiff < 0)] <- 0

ggplot(lmap, aes(GenomicPosition, cMDiff)) + geom_point()

write.table(lmap, 
            paste("../../results/2_Merged_map_", AnalysisSuffix, ".txt", sep = ""), 
            sep = "\t", quote = F, row.names = F)

#~~ Create a table of the maximum values of each chromosome

maxvals <- data.frame(maxcM = tapply(lmap$new.cM, lmap$Chr, max),
                      maxcM.male = tapply(lmap$cMPosition.Male, lmap$Chr, max),
                      maxcM.female = tapply(lmap$cMPosition.Female, lmap$Chr, max),
                      maxGenome = tapply(lmap$GenomicPosition, lmap$Chr, max),
                      Chr = unique(lmap$Chr),
                      nsnps = tapply(lmap$Order, lmap$Chr, max))


#~~ Examine the correspondence between linkage and physical map lengths

ggplot(maxvals, aes(maxGenome/1000000, maxcM)) +
  stat_smooth(method = "lm", alpha = 0.2) +
  geom_text(aes(label = Chr), size = 5, fontface = "bold") +
  defopts +
  scale_x_continuous(breaks = c(seq(50, 300, 50))) +
  labs(x = "Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(maxvals, aes(maxGenome/1000000, maxcM.male)) +
  stat_smooth(method = "lm", alpha = 0.2) +
  geom_text(aes(label = Chr), size = 5, fontface = "bold") +
  defopts +
  scale_x_continuous(breaks = c(seq(50, 300, 50))) +
  labs(x = "Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(maxvals, aes(maxGenome/1000000, maxcM.female)) +
  stat_smooth(method = "lm", alpha = 0.2) +
  geom_text(aes(label = Chr), size = 5, fontface = "bold") +
  defopts +
  scale_x_continuous(breaks = c(seq(50, 300, 50))) +
  labs(x = "Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(subset(maxvals, Chr != 27), aes(as.factor(Chr), (maxcM.male - maxcM.female)/maxGenome)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "Difference in Linkage Map Lengths Corrected for Length (cM)")


#~~ per chromosome...

ggplot(lmap, aes(GenomicPosition, new.cM)) +
  geom_point(alpha = 0.2) +
  defopts +
  labs(x = "Genomic Position (MB)", y = "Linkage Map Position (cM)") +
  facet_wrap(~Chr, scales = "free")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Determine the Order of first and last informative SNPs on each chromosome  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Determine the Informative Start and Stop Positions

rectab$First.Inf.Order <- unlist(lapply(rectab$data, InfLengthFunc))
rectab$Last.Inf.Order <-  unlist(lapply(rectab$data, function(x) InfLengthFunc(x, position = "Last")))
rectab$First.Inf.Order.Cons <- unlist(lapply(rectab$data, InfLengthFunc.Cons))
rectab$Last.Inf.Order.Cons <-  unlist(lapply(rectab$data, function(x) InfLengthFunc.Cons(x, position = "Last")))

rectab$RecombCount.Cons <- unlist(lapply(rectab$data, ConservativeRecCountFunc))


beep()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Determine the informative length of the chromosome using map information   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ merge with map information

rectab <- join(rectab, maxvals[,c(4:6)])

names(rectab)[which(names(rectab) %in% c("maxGenome", "nsnps"))] <- c("Chromosome.Length", "Chromosome.SNP.Count")


lmap$First.Inf.Order <- lmap$Order
lmap$Last.Inf.Order  <- lmap$Order

maptabstart <- subset(lmap, select = c(Chr, GenomicPosition, First.Inf.Order))
maptabstop  <- subset(lmap, select = c(Chr, GenomicPosition, Last.Inf.Order ))

names(maptabstart)[2] <- "First.Inf.Pos"
names(maptabstop)[2] <- "Last.Inf.Pos"

rectab <- join(rectab, maptabstart)
rectab <- join(rectab, maptabstop)

#~~ Determine the proportion of the chromosome captured between the first and last markers
rectab$Prop.Inf.Chromosome <- (rectab$Last.Inf.Pos - rectab$First.Inf.Pos)/rectab$Chromosome.Length

rectab$Prop.Inf.SNPs <- rectab$No.Inf.Loci/rectab$Chromosome.SNP.Count

#~~ Determine the recombination Rate per chromosome based on the length of the informative chromosome.

rectab$Recomb.Rate <- rectab$RecombCount * 1e6/(rectab$Last.Inf.Pos - rectab$First.Inf.Pos)

rm(maptabstart, maptabstop, i)

#~~ Write to file

write.table(maxvals, paste("../../results/2_MaxVals_", AnalysisSuffix, ".txt", sep = ""), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Determine the origin of informative SNPs and count them                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Not required - see superceded versions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Extract information about the individual in which recombination took place #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(crimapped)
names(crimapped)[1] <- "id"

rectab <- join(rectab, crimapped)

#~~ define which individual the Recombination actulaly took place in



rectab$Offspring.ID <- unlist(lapply(rectab$Family, function(x) strsplit(x, split = "_")[[1]][3]))


rectab <- subset(rectab, Offspring.ID == id)
rectab <- subset(rectab, No.Inf.Loci > 0)

head(rectab)

rectab$RRID <- NA
rectab$RRID[which(rectab$parent == "MOTHER")] <- rectab$MOTHER[which(rectab$parent == "MOTHER")]
rectab$RRID[which(rectab$parent == "FATHER")] <- rectab$FATHER[which(rectab$parent == "FATHER")]

#~~ What is the recombination rate?

for(i in 1:ncol(rectab)) attributes(rectab[,i]) <- NULL

rectab$Inf.Chr.Length <- rectab$Last.Inf.Pos - rectab$First.Inf.Pos
rectab$RR <- rectab$RecombCount/rectab$Inf.Chr.Length


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Add information about the ID, it's offspring, it's age, population data, etc. #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(rectab)

head(basedata)

#~~ merge Offspring and RRID sex and birth year information

basedata <- subset(basedata, select = c(ID, SEX, BIRTHYEAR))

names(basedata) <- c("RRID", "RRID.SEX", "RRID.BYEAR")
rectab <- join(rectab, basedata)

names(basedata) <- c("Offspring.ID", "Offspring.SEX", "Offspring.BYEAR")
rectab <- join(rectab, basedata)

#~~ determine CAPAGE of RRID

rectab$RRID.CAPAGE <- rectab$Offspring.BYEAR - rectab$RRID.BYEAR


write.table(rectab,
            paste("../../results/2_IndivRRRaw_", AnalysisSuffix, ".txt", sep = ""), 
            row.names = F, sep = "\t", quote = F)
