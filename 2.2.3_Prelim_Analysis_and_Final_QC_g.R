
# QC and Data Preparation
# Author: Susan Johnston



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ load functions and libraries

library(beepr)
library(ggplot2)
library(reshape)
library(lme4)
library(GenABEL)
library(plyr)


source("r/recoderFunc.R")
source("r/makeGRM.R")


options(stringsAsFactors = F)

defopts <- theme(axis.text.x  = element_text (size = 14, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 14, vjust = 0.7),
                 axis.title.y = element_text (size = 14, angle = 90, vjust = 0.2),
                 axis.title.x = element_text (size = 14, vjust = 0.2),
                 strip.background = element_blank(),
                 legend.position = "top")

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in genomic map positions

lmap <- read.table(paste0("results/2_Merged_map_", AnalysisSuffix, ".txt"), header = T)

#~~ max values per chromosome

maxvals <- read.table(paste0("results/2_MaxVals_", AnalysisSuffix, ".txt"), header = T)

#~~ read in recombination information

rectab <- read.table(paste0("results/2_IndivRRDoubClean_", AnalysisSuffix,".txt"), header = T, stringsAsFactors = F)
recsumm <- read.table(paste0("results/2_TotalSummIndivRRDoubClean_", AnalysisSuffix,".txt"), header = T, stringsAsFactors = F)

#~~ read in pedigrees

pedigree <- read.table("data/pedigree_20130920.txt", header = T)
crimapped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T)

#~~ read in genomic relatedness information

memory.limit(size = 800000000) #ESSENTIAL OTHERWISE YOU RUN INTO MEMORY ISSUES

grm.gcta <- read.table("gcta/150129_autoGRM_adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.gcta <- read.table("gcta/150129_autoGRM_adj.grm.id")  # CONTAINS ID LIST

head(grm.gcta)
head(ids.gcta)

str(grm.gcta)
str(ids.gcta)

#~~ subset the GRM to include only individuals within crimapped

ids.gcta <- ids.gcta[which(ids.gcta$V2 %in% c(crimapped$ANIMAL, crimapped$FATHER, crimapped$MOTHER)),]
grm.gcta <- grm.gcta[which(grm.gcta$V1 %in% ids.gcta$V1 & grm.gcta$V2 %in% ids.gcta$V1),-3]

#~~ Read in SNP data

sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")

summary(qtscore(NewHorn2, data = sheepabel))


#~~ read in inbreeding information

inbreedingtab <- read.table("gcta/20150129_autosomal_IBD.ibc", header = T)

inbreedingtab <- inbreedingtab[,c("IID", "Fhat3")]
names(inbreedingtab) <- c("RRID", "RRID.Fhat3")

rectab  <- join(rectab , inbreedingtab)
recsumm <- join(recsumm, inbreedingtab)

names(inbreedingtab) <- c("Offspring.ID", "Offspring.Fhat3")

rectab  <- join(rectab , inbreedingtab)
recsumm <- join(recsumm, inbreedingtab)

#~~ transform RecombRate to be in terms of MB rather than bases

rectab$RecombRate           <- rectab$RecombRate * 1e6
recsumm$MeanRRincNonRecombs <- recsumm$MeanRRincNonRecombs * 1e6

#~~ give Recombination Rate in recsumm in terms of informative genome

recsumm$TotalRecombRate <- recsumm$TotalRecombCount/recsumm$TotalInfChrLenIncNonRecombs
recsumm$TotalRecombRate <- recsumm$TotalRecombRate * 1e6


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### 1. Chromosome and mapping information      #                             
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

maxvals$maxGenome2 <- maxvals$maxGenome/1000000

#~~ Physical vs. linkage map length
ggplot(maxvals, aes(maxGenome2, maxcM)) +
  stat_smooth(method = "lm", alpha = 0.2) +                     # formula = y ~ log(x)
  geom_text(aes(label = Chr), size = 5, fontface = "bold") +
  defopts +
  scale_x_continuous(breaks = c(seq(50, 300, 50))) +
  labs(x = "Chromosome Length (MB)",
       y = "Linkage Map Length (cM)",
       title = paste0("Adj.R.Squared = ", summary(lm(maxcM ~ maxGenome2, data = maxvals))$adj.r.squared))

#~~ Chromosome by chromosome
ggplot(lmap, aes(GenomicPosition/1000000, cMPosition)) +
  geom_point(size = 2, alpha = 0.1) +
  defopts +
  facet_wrap(~ Chr, scales = "free") +
  labs(x = "Genomic Position (MB)", y = "Linkage Map Position (cM)")

#~~ Chromosome by chromosome by sex

meltmap <- melt(lmap[,c("Chr", "cMPosition.Female", "cMPosition.Male", "GenomicPosition")], id.vars = c("Chr", "GenomicPosition"))
head(meltmap)
meltmap$variable <- gsub("cMPosition.", "", meltmap$variable)
  
ggplot(meltmap, aes(GenomicPosition/1000000, value, colour = variable)) +
  geom_point(size = 2, alpha = 0.1) +
  defopts +
  facet_wrap(~ Chr, scales = "free") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Genomic Position (MB)", y = "Linkage Map Position (cM)", col = "")

#~~ Recombination fraction difference

head(lmap)

lmap$rDiff <- lmap$Male.r - lmap$Female.r

ggplot(lmap, aes(GenomicPosition, rDiff)) +
  geom_line(alpha = 0.3) +
  defopts +
  facet_wrap(~ Chr, scales = "free_x") +
  labs(x = "Genomic Position (MB)", y = "Difference in r (Male - Female)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### 2. Recombination Information               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Basic data summaries using data across all chromosomes (i.e. not genome-wide summary data)

ggplot(rectab, aes(RecombCount.v2))      + geom_histogram(binwidth = 1   , col = "grey") + defopts + labs(x = "Number of Crossovers")
ggplot(rectab, aes(Prop.Inf.SNPs))       + geom_histogram(binwidth = 0.01, col = "grey") + defopts + labs(x = "Proportion of Informative SNPs")
ggplot(rectab, aes(Prop.Inf.Chromosome)) + geom_histogram(binwidth = 0.01, col = "grey") + defopts + labs(x = "Proportion of Chr with Informative SNPs")
ggplot(rectab, aes(RecombRate))          + geom_histogram(                 col = "grey") + defopts + labs(x = "Recombination Rate (per MB)")

ggplot(rectab, aes(Prop.Inf.SNPs)) + 
  geom_histogram(binwidth = 0.01, col = "grey") + 
  defopts + 
  labs(x = "Proportion of Informative SNPs") +
  facet_wrap(~ Chr)


#~~ by chromosome

ggplot(rectab, aes(factor(Chr), RecombCount.v2))    + geom_boxplot() + defopts + labs(x = "Chromosome", y = "Number of Crossovers")
ggplot(rectab, aes(factor(Chr), RecombRate))        + geom_boxplot() + defopts + labs(x = "Chromosome", y = "Recombination Rate (Xover per MB)")


#~~ Non recombinant chromosomes

nonrectab <- data.frame(table(rectab$RecombCount.v2, rectab$Chr))
head(nonrectab)
names(nonrectab) <- c("RecombCount", "Chr", "Freq")

nonrectab$Freq <- nonrectab$Freq/nrow(recsumm)
nonrectab$NonRecomb <- nonrectab$RecombCount
nonrectab$NonRecomb[which(as.numeric(as.character(nonrectab$RecombCount)) > 0)] <- 1

nonrectab <- droplevels(nonrectab)

ggplot(nonrectab, aes(Chr, Freq, fill = NonRecomb)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  defopts


# Basic data summaries using data summarised across all chromosomes in a single
# reproductive event. These values include non-recombinant chromosomes.

ggplot(recsumm, aes(TotalRecombCount))                    + geom_histogram(binwidth = 1,     col = "grey") + defopts + labs(x = "Total Number of Crossovers")
ggplot(recsumm, aes(TotalRecombCount))                    + geom_histogram(binwidth = 1,     col = "grey") + defopts + labs(x = "Total Number of Crossovers") + facet_wrap(~RRID.SEX)
ggplot(recsumm, aes(TotalInfLoci))                        + geom_histogram(binwidth = 100,   col = "grey") + defopts + labs(x = "Total Number of Informative Loci")
ggplot(recsumm, aes(MeanPropInfLoci))                     + geom_histogram(binwidth = 0.005, col = "grey") + defopts + labs(x = "Mean Proportion of Informative Loci")
ggplot(recsumm, aes(MeanPropChr))                         + geom_histogram(binwidth = 0.005, col = "grey") + defopts + labs(x = "Mean Proportion of the Chromosome Represented")
ggplot(recsumm, aes(UninfChrCount))                       + geom_histogram(                  col = "grey") + defopts + labs(x = "Number of Uninformative Chromosomes")
ggplot(recsumm, aes(TotalInfChrLenIncNonRecombs/1000000)) + geom_histogram(binwidth = 10,    col = "grey") + defopts + labs(x = "Total Informative Chromosome Length (MB)")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Quality control on recombination measures          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ What is the total possible informative length of the genome?

informtab <- data.frame(MinPos = tapply(lmap$GenomicPosition, lmap$Chr, min),
                        MaxPos = tapply(lmap$GenomicPosition, lmap$Chr, max))

informtab$ChrTotal <- informtab$MaxPos - informtab$MinPos

recsumm$TotalPropInfGenome <- recsumm$TotalInfChrLenIncNonRecombs/sum(as.numeric(informtab$ChrTotal))

ggplot(recsumm, aes(TotalPropInfGenome)) + 
  geom_histogram(binwidth = 0.001,    col = "grey") + 
  defopts + 
  labs(x = "Proportion of Informative Genome")

ggplot(rectab, aes(Prop.Inf.Chromosome)) + 
  geom_histogram(binwidth = 0.01, col = "grey") + 
  defopts + 
  labs(x = "Proportion of Chr with Informative SNPs")

table(table(rectab[which(rectab$Prop.Inf.Chromosome < 0.8),"Family"]))

#~~ Are there any really inbred individuals?

ggplot(recsumm, aes(RRID.Fhat3)) + geom_histogram(col = "grey") + defopts
ggplot(recsumm, aes(RRID.Fhat3, TotalRecombRate)) + geom_point() + defopts + stat_smooth(method = "lm")
ggplot(recsumm, aes(RRID.Fhat3, TotalInfLoci)) + geom_point() + defopts
ggplot(recsumm, aes(RRID.Fhat3, TotalPropInfGenome)) + geom_point() + defopts


#~~ Remove IDs with less than 0.9 of informative genome

removefam <- recsumm[which(recsumm$TotalPropInfGenome < 0.9), "Family"]

#~~ remove these families from analysis datasets

recsumm <- subset(recsumm, !Family %in% removefam)
rectab <- subset(rectab, !Family %in% removefam)

#~~ Remove IDs with less than 8,000 informative SNPs

removefam <- recsumm[which(recsumm$TotalInfLoci < 7000), "Family"]

#~~ remove these families from analysis datasets

recsumm <- subset(recsumm, !Family %in% removefam)
rectab <- subset(rectab, !Family %in% removefam)

#~~ Remove Ids where an individual occurs twice in their nuclear family

famedit <- data.frame(table(crimapped$Family, crimapped$ANIMAL))
head(famedit)
famedit <- subset(famedit, Freq > 1)

#~~ remove these families from analysis datasets

rectab <- subset(rectab, !Family %in% famedit$Var1)
recsumm <- subset(recsumm, !Family %in% famedit$Var1)


#~~ Clean up workspace

rm(informtab, removefam, famedit)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### 4. Determine relatedness between individuals.....   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

###~~ MATE ID AND INBREEDING ~~###

#~~ Create a table with RRID and MateID for each OffspringID

temped <- subset(pedigree, ANIMAL %in% recsumm$Offspring.ID)

names(temped) <- c("Offspring.ID", "RRID", "MateID")

temped2 <- temped[,c(1, 3, 2)]

head(temped2)

names(temped2) <- c("Offspring.ID", "RRID", "MateID")

matetab <- rbind(temped, temped2)

#~~ Merge with rectab and recsumm

recsumm <- join(recsumm, matetab)
rectab <- join(rectab, matetab)

#~~ Clean up workspace

rm(temped, temped2, matetab)

#~~ Add statistics on the mate

head(inbreedingtab)
names(inbreedingtab) <- c("MateID", "Mate.Fhat3")

recsumm <- join(recsumm, inbreedingtab)
rectab  <- join(rectab,  inbreedingtab)


####~~ RELATEDNESS OF FAMILY UNIT ~~####

head(grm.gcta)
head(ids.gcta)

famtab <- data.frame(Family = unique(recsumm$Family),
                     Mean.FamR.wids = NA,
                     Mean.FamR.woids = NA)

system.time({
  for(i in 1:nrow(famtab)){
    
    if(i %in% seq(1, nrow(famtab), 10)) print(paste("Analysing Family", i, "of", nrow(famtab)))
    
    fam <- famtab$Family[i]
    pedids <- crimapped[which(crimapped$Family == fam),"ANIMAL"]
    gctaids <- ids.gcta[which(ids.gcta$V2 %in% pedids),"V1"]
    
    famgcta <- grm.gcta[which(grm.gcta$V1 %in% gctaids & grm.gcta$V2 %in% gctaids),]
    
    famtab$Mean.FamR.wids[i]  <- mean(famgcta$V4)
    famtab$Mean.FamR.woids[i] <- mean(famgcta[which(famgcta$V1 != famgcta$V2), "V4"])
    
    rm(fam, pedids, gctaids, famgcta)
  }
})

write.table(famtab, paste0("results/2_Family_Genomic_Relatedness_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)

famtab <- read.table(paste0("results/2_Family_Genomic_Relatedness_", AnalysisSuffix, ".txt"), header = T)

ggplot(famtab, aes(Mean.FamR.wids)) + geom_histogram(binwidth = 0.005, col = "grey")

ggplot(famtab, aes(Mean.FamR.woids)) + geom_histogram(binwidth = 0.005, col = "grey")

recsumm <- join(recsumm, famtab)
rectab  <- join(rectab,  famtab)

ggplot(recsumm, aes(TotalInfLoci, Mean.FamR.wids)) + geom_point()
ggplot(recsumm, aes(TotalInfLoci, Mean.FamR.woids)) + geom_point() + stat_smooth(method = "lm")


#~~ fix crimapped
crimapped <- crimapped[which(crimapped$Family %in% recsumm$Family),]

write.table(crimapped, paste0("results/2_FamilyPedigree_FullClean_", AnalysisSuffix, ".txt"),    row.names = F, sep = "\t", quote = F)
write.table(rectab,    paste0("results/2_IndivRR_FullClean_", AnalysisSuffix, ".txt"),    row.names = F, sep = "\t", quote = F)
write.table(recsumm,   paste0("results/2_TotalSummIndivRR_FullClean_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)



##########################################################################################################
######~~~~~~~~~~~~~~~~~~~ AFTER THIS POINT IS PURELY EXPLORATOR ANALYSES - NOW SUPERCEDED ~~~~~~~~~~~~~~#
##########################################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### 5. What else could affect the TotalInfLoci measure? #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

famorder <- readLines("crimap/crimap_f/chr24f.gen")
famorder <- famorder[grep("Offspring_", famorder)]

famorder <- data.frame(Family = famorder)
famorder$FamOrder <- 1:nrow(famorder)

recsumm <- join(recsumm, famorder)

ggplot(recsumm, aes(FamOrder, TotalInfLoci, col = RRID.SEX)) + geom_point() + defopts + scale_colour_brewer(palette = "Set1") + labs(x = "Total Number of Crossovers", y = "Total Number of Informative Loci")
ggplot(recsumm, aes(TotalRecombCount, TotalInfLoci, col = RRID.SEX)) + geom_point() + defopts + scale_colour_brewer(palette = "Set1") + labs(x = "Total Number of Crossovers", y = "Total Number of Informative Loci")
ggplot(recsumm, aes(TotalRecombRate,  TotalInfLoci, col = RRID.SEX)) + geom_point() + defopts + scale_colour_brewer(palette = "Set1") + labs(x = "Recombination Rate (Per MB)", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(as.factor(RRID.BYEAR), TotalInfLoci)) + 
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  labs(x = "RRID BIRTH YEAR", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(as.factor(Offspring.BYEAR), TotalInfLoci)) +
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  labs(x = "Offspring BIRTH YEAR", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(RRID.BYEAR, TotalInfLoci)) + 
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_point(notch = T) + 
  defopts + 
  stat_smooth(method = "lm") +
  labs(x = "RRID BIRTH YEAR", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(Offspring.BYEAR, TotalInfLoci)) +
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_point(notch = T) + 
  defopts + 
  stat_smooth(method = "lm") +
  labs(x = "Offspring BIRTH YEAR", y = "Total Number of Informative Loci")




ggplot(recsumm, aes(NonRecombChrCount)) + geom_histogram(binwidth = 1,     col = "grey") + defopts + labs(x = "Number of Non-Recombinant Chromosomes")
mean(recsumm$NonRecombChrCount)


ggplot(recsumm, aes(RRID.SEX, TotalRecombCount)) + geom_boxplot(notch = T) + defopts + labs(x = "Sex of RRID", y = "Total Number of Crossovers")

ggplot(recsumm, aes(RRID.SEX, MeanRRincNonRecombs)) + 
  geom_boxplot(notch = T, width = 0.4) + 
  defopts + 
  labs(x = "Sex", y = "Genome-wide Recombination Rate (Xovers per MB)") 

str(glm(MeanRRincNonRecombs ~ RRID.SEX, data = recsumm))

ggplot(recsumm, aes(Offspring.SEX, TotalRecombCount)) + geom_boxplot(notch = T) + defopts + labs(x = "Sex of Offspring", y = "Total Number of Crossovers")

ggplot(recsumm, aes(factor(RRID.CAPAGE), TotalRecombCount, fill = RRID.SEX)) +
  geom_boxplot(notch = T) + 
  defopts + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Age of RRID", y = "Total Number of Crossovers")+
  facet_wrap(~RRID.SEX, scales = "free_x")

recsumm$CAPAGE2 <- recsumm$RRID.CAPAGE
recsumm$CAPAGE2[which(recsumm$RRID.SEX == "Female")] <- recsumm$CAPAGE[which(recsumm$RRID.SEX == "Female")] + 0.2


ggplot(recsumm, aes(CAPAGE2, MeanRRincNonRecombs, colour = RRID.SEX)) +
  geom_point(alpha = 0.5) + 
  stat_smooth(method = "lm") +
  defopts + 
  scale_x_continuous(breaks = seq(1, 15, 1)) +
  scale_colour_brewer(palette = "Set1") + 
  labs(x = "Age at Gamete Tranfer", y = "Genome-wide RR (Xovers per MB)")

summary(glm(MeanRRincNonRecombs ~ RRID.SEX * RRID.CAPAGE, data = recsumm))


ggplot(recsumm, aes(as.factor(RRID.BYEAR), TotalRecombCount, fill = RRID.SEX)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  labs(x = "RRID Birth Year", y = "Total Number of Crossovers") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~RRID.SEX, scales = "free_x")


ggplot(recsumm, aes(as.factor(Offspring.BYEAR), TotalRecombCount, fill = RRID.SEX)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  labs(x = "Offspring Birth Year", y = "Total Number of Crossovers")+
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~RRID.SEX, scales = "free_x")

ggplot(recsumm, aes(as.factor(RRID.BYEAR), TotalInfLoci)) + 
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  labs(x = "RRID BIRTH YEAR", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(as.factor(Offspring.BYEAR), TotalInfLoci)) +
  theme(axis.text.x  = element_text (size = 16, vjust = 0, angle = 270)) + 
  geom_boxplot(notch = T) + 
  defopts + 
  labs(x = "Offspring BIRTH YEAR", y = "Total Number of Informative Loci")

source("C:/Users/Susan Johnston/Desktop/R Functions/countIF.R")
recsumm$RRID.Count <- countIF(recsumm$RRID)


table(recsumm$RRID.Count)/as.numeric(names(table(recsumm$RRID.Count)))

ggplot(recsumm, aes(RRID.Count, TotalInfLoci)) + geom_point() + defopts + labs(x = "RRID Count", y = "Total Number of Informative Loci")

ggplot(recsumm, aes(RRID.Count, TotalRecombCount, col = RRID.SEX)) + 
  geom_point(alpha = 0.3) + 
  stat_smooth(method = "lm") + 
  scale_colour_brewer(palette = "Set1") +
  defopts +
  facet_wrap(~RRID.SEX, scale = "free_x") +
  labs(x = "RRID Count", y = "Total Number of Crossovers")




summary(lmer(TotalRecombCount ~ RRID.Count * RRID.SEX + (1|RRID), data = recsumm))
summary(lmer(TotalRecombCount ~ RRID.CAPAGE * RRID.SEX + (1|RRID), data = recsumm))




library(asreml)

head(pedigree)
names(pedigree)[1] <- "RRID"
pedigree$RRID <- as.factor(pedigree$RRID)
pedigree$MOTHER <- as.factor(pedigree$MOTHER)
pedigree$FATHER <- as.factor(pedigree$FATHER)

recsumm <- merge(recsumm, pedigree, all.x = T)

recsumm$RRID <- as.factor(recsumm$RRID)
recsumm$MOTHER <- as.factor(recsumm$MOTHER)
recsumm$FATHER <- as.factor(recsumm$FATHER)


ainv <- asreml.Ainverse(pedigree)$ginv
str(recsumm)

model1 <- asreml(fixed = TotalRecombCount ~ factor(RRID.SEX),
                 random = ~ped(RRID) + ide(RRID) + RRID.BYEAR + Offspring.BYEAR + ide(MOTHER),
                 data = recsumm,
                 ginverse =  list(RRID = ainv, MOTHER = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6)


model1.fem <- asreml(fixed = TotalRecombCount ~ factor(RRID.SEX),
                     random = ~ped(RRID) + ide(RRID) + RRID.BYEAR + Offspring.BYEAR + ide(MOTHER),
                     data = subset(recsumm, RRID.SEX == "Female"),
                     ginverse =  list(RRID = ainv, MOTHER = ainv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6)


model1.mal <- asreml(fixed = TotalRecombCount ~ factor(RRID.SEX),
                     random = ~ped(RRID) + ide(RRID) + RRID.BYEAR + Offspring.BYEAR + ide(MOTHER),
                     data = subset(recsumm, RRID.SEX == "Male"),
                     ginverse =  list(RRID = ainv, MOTHER = ainv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6)


x <- summary(model1, all = T)

source("C:/Users/Susan Johnston/Desktop/R Functions/ASReml.EstEffects.R")

model2 <- asreml(fixed = MeanRRincNonRecombs ~ factor(RRID.SEX),
                 random = ~ped(RRID) + ide(RRID) + RRID.BYEAR + Offspring.BYEAR,
                 data = recsumm,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6)

x.2 <- summary(model2, all = T)

memory.limit(size = 800000000) #ESSENTIAL OTHERWISE YOU RUN INTO MEMORY ISSUES

grm.gcta <- read.table("gcta/allnewgrmadj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.gcta <- read.table("gcta/allnewgrmadj.grm.id")  # CONTAINS ID LIST
# # 

source("C:/Users/Susan Johnston/Desktop/R Functions/makeGRM.R")

recsumm$ID <- factor(recsumm$RRID)
recsumm$ID2 <- factor(recsumm$RRID)


x.all <- makeGRM(grm.gcta,ids.gcta,recsumm, which(names(recsumm) == "ID"))
recsumm2 <- droplevels(subset(recsumm, ID %in% dimnames(x.all$grm)[[1]]))

lambda(model1.2)


table(recsumm$ID %in% dimnames(x.all$grm)[[1]])

recsumm2$CAPAGE <- as.numeric(recsumm2$RRID.CAPAGE)
recsumm2$RRID.BYEAR <- as.factor(recsumm2$RRID.BYEAR)
recsumm2$RRID.SEX <- as.factor(recsumm2$RRID.SEX)

recsumm2$Offspring.BYEAR <- as.factor(recsumm2$Offspring.BYEAR)

names(recsumm2)


str(recsumm2)


grm.model <- asreml(fixed  =  TotalRecombCount ~ CAPAGE + RRID.SEX + RRID.Fhat3 + TotalInfLoci,
                    random   =  ~ giv(ID) + ide(ID) + RRID.BYEAR + Offspring.BYEAR, 
                    data     =  recsumm2,
                    ginverse =  list(ID = x.all$grm),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)


grm.model.M <- asreml(fixed  =  MeanRRincNonRecombs ~ CAPAGE + RRID.Fhat3,
                      random   =  ~ giv(ID) + ide(ID) + RRID.BYEAR + Offspring.BYEAR, 
                      data     =  subset(recsumm2, RRID.SEX == "Male"),
                      ginverse =  list(ID = x.all$grm),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

grm.model.F <- asreml(fixed  =  MeanRRincNonRecombs ~ CAPAGE + RRID.Fhat3,
                      random   =  ~ giv(ID) + ide(ID) + RRID.BYEAR + Offspring.BYEAR, 
                      data     =  subset(recsumm2, RRID.SEX == "Female"),
                      ginverse =  list(ID = x.all$grm),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

str(ids.gcta)


grm.modelextra <- asreml(fixed  =  MeanRRincNonRecombs ~ CAPAGE + RRID.SEX + RRID.Fhat3 + TotalInfLoci,
                         random   =  ~ ide(ID) + RRID.BYEAR + Offspring.BYEAR, 
                         data     =  recsumm2,
                         ginverse =  list(ID = x.all$grm),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)



x$coef.fixed
ASReml.EstEffects(model1)

x <- summary(grm.model, all = T)
x.M <- summary(grm.model.M, all = T)
x.F <- summary(grm.model.F, all = T)

x$coef.fixed
x.M$coef.fixed
x.F$coef.fixed

asreml::wald.asreml(grm.model)
asreml::wald.asreml(grm.model.M)
asreml::wald.asreml(grm.model.F)


ASReml.EstEffects(grm.model)
ASReml.EstEffects(grm.model.M)
ASReml.EstEffects(grm.model.F)


ASReml.EstEffects(model1.fem)
ASReml.EstEffects(model1.mal)


recsumm3 <- subset(recsumm, TotalInfLoci < 9000)
ggplot(recsumm3, aes(TotalInfLoci, RRID.Fhat3)) + 
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  defopts +
  labs(x = "Total Number of Informative Loci", y = "Genomic ^F")

ggplot(recsumm3, aes(TotalInfLoci, MeanRRincNonRecombs)) + 
  geom_point(alpha = 0.8) +
  defopts +
  stat_smooth(method = "lm") +
  labs(x = "Total Number of Informative Loci", y = "Genome-wide Recombination Rate (per MB)")





### GWAS of Recombination Count.


weetab <- data.frame(cbind(as.numeric(grm.model$residuals), as.character(recsumm2[,"RRID"])))
names(weetab) <- c("Residual2", "RRID")

weetab$Residual2 <- as.numeric(weetab$Residual2)


head(weetab)

weetab <- data.frame(Residual2 = tapply(weetab$Residual2, as.factor(weetab$RRID), mean))

head(weetab)
weetab$id <- row.names(weetab)
str(model1)


rridres <- data.frame(grm.model$coefficients$random)
rridres$Info <- row.names(rridres)
rridres <- rridres[grep("giv(ID)_", rridres$Info, fixed = T),]

rridres$Info <- gsub("giv(ID)_", "", as.character(rridres$Info), fixed = T)

names(rridres) <- c("RRIDGBLUP", "id")

rridres <- rridres[which(rridres$id %in% recsumm$RRID),]

sheepabel2 <- load.gwaa.data(phe = "data/1_GenABEL_Basedata20140630.txt",
                             gen = "data/1_GenABEL_sheepabel_pedQC20140630.gen")

summary(qtscore(MeanResid.Leng, data = sheepabel2))


sheepabel2 <- add.phdata(sheepabel2, rridres)
sheepabel2 <- add.phdata(sheepabel2, weetab)



model1.1 <- qtscore(RRIDGBLUP ~ sex, data = sheepabel2)
model1.2 <- qtscore(Residual2 ~ sex, data = sheepabel2)

names(phdata(sheepabel2))


source("C:/Users/Susan Johnston/Desktop/R Functions/GenABELPlotFunctions.R")
source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")

bonf = 0.05/nsnps(sheepabel2)
FullSummary(model1.1, bonf)
FullSummary(model1.2, bonf)


rectab.chr3 <- subset(rectab, Chr == 3)
sheepabel2.chr3 <- sheepabel2[,which(chromosome(sheepabel2) == 3)]
rectab.chr3 <- subset(rectab.chr3, select = c(Recomb.Rate, RRID))
names(rectab.chr3)[2] <- "id"
rectab.chr3 <- data.frame(tapply(rectab.chr3$Recomb.Rate, rectab.chr3$id, mean))

head(rectab.chr3)
names(rectab.chr3) <- "Chr3.RR"
rectab.chr3$id <- row.names(rectab.chr3)

sheepabel2.chr3 <- add.phdata(sheepabel2.chr3, rectab.chr3)

test <- qtscore(Chr3.RR ~ sex, sheepabel2.chr3)
summary(test)

test2 <- results(test)
head(test2)

test2 <-   test2[with(test2, order(Position)), ]

test2 <- data.frame(test2)
test2$Diff <-   c(0,diff(test2$Position))
test2$Diff2 <-   c(diff(test2$Position), 0)
test2$MeanDiff <- rowSums(test2[,c("Diff", "Diff2")])


head(test2)


plot(test2$MeanDiff, -log10(test2$Pc1df))

head(rectab)





#~~ 



peridsumm <- data.frame(IDCount     = tapply(recsumm$RRID, recsumm$RRID, length),
                        MeanInfLoci = tapply(recsumm$TotalInfLoci, recsumm$RRID, mean))

hist(peridsumm$MeanInfLoci)
hist(peridsumm$IDCount)

plot(peridsumm$MeanInfLoci, peridsumm$IDCount)


recsumm$UniqueID2 <- paste(recsumm$RRID, recsumm$Offspring.ID, sep = "_")

rectab$UniqueID2 <- paste(rectab$RRID, rectab$Offspring.ID, sep = "_")

highinf <- recsumm[which(recsumm$TotalInfLoci > 11000),"UniqueID2"]
highinf.rectab <- subset(rectab, UniqueID2 %in% highinf)


ggplot(recsumm, aes(TotalRecombCount, TotalInfLoci)) + geom_point()
