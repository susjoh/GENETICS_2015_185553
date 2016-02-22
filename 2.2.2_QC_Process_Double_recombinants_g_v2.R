# Conduct Quality Control and Process the double recombinants
# Author: Susan Johnston


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#~~ load functions and libraries

library(beepr)
library(ggplot2)
library(plyr)

source("r/recoderFunc.R")

save.test <- function(data) write.table(data, "test.txt", row.names = F, quote = F, sep = "\t")

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in data

rectab   <- read.table(paste0("results/2_IndivRRRaw_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)
lmap     <- read.table(paste0("results/2_Merged_map_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)
pedigree <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)
pseudoautosomal <- read.table("results/1_Pseudoautosomal_SNPs_in_X.txt", header = T, stringsAsFactors = F)

#~~ remove families with duplicate IDs

famedit <- data.frame(table(pedigree$Family, pedigree$ANIMAL))
head(famedit)

famedit <- subset(famedit, Freq > 1)

pedigree <- subset(pedigree, !Family %in% famedit$Var1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Are there some individuals with particularly unusual recombination rate?   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ create a vector of RRID combined with the offspring ID.

rectab$RRID.Offspring.ID <- paste(rectab$RRID, rectab$Offspring.ID, sep = "_")

table(table(rectab$RRID.Offspring.ID))

#~~ Create a table of Total Recombination Count for each gamete transfer

recomb.summ <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$RRID.Offspring.ID, sum),
                          TotalRecombCount.Cons = tapply(rectab$RecombCount.Cons, rectab$RRID.Offspring.ID, sum))
recomb.summ$RRID.Offspring.ID <- row.names(recomb.summ)
recomb.summ$RRID <- unlist(lapply(recomb.summ$RRID.Offspring.ID, function(string) unlist(strsplit(string, split = "_")[[1]][1])))
recomb.summ$Offspring.ID <- unlist(lapply(recomb.summ$RRID.Offspring.ID, function(string) unlist(strsplit(string, split = "_")[[1]][2])))

head(recomb.summ)
hist(recomb.summ$TotalRecombCount)
table(recomb.summ$TotalRecombCount)

ggplot(recomb.summ, aes(TotalRecombCount, TotalRecombCount.Cons)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm")


#~~ Do some chromatids have an unusually high number of crossovers?

rectab$Chr <- factor(rectab$Chr)

ggplot(rectab, aes(Chr, RecombCount)) + geom_boxplot()

ggplot(rectab, aes(RecombCount)) + geom_histogram(binwidth = 1) + facet_wrap(~Chr, scales = "free")

rectab$data[which(rectab$RecombCount > 7 & !rectab$Chr %in% c(1, 2, 3))]

tapply(rectab$RecombCount, rectab$Chr, mean)
tapply(rectab$RecombCount, rectab$Chr, sd)*5

#~~ remove unreasonable recombination events

rectab <- subset(rectab, !Family %in%  rectab[which(rectab$RecombCount > 7 & !rectab$Chr %in% c(1, 2, 3)),"Family"])

#~~ identify rectab entries with singletons

SingletonStatusFunc <- function(vector){
  
  unlist(lapply(vector, function(string){
    
    y <- unlist(strsplit(string, split = ""))
    y <- y[which(y != "-")]
    y[which(y %in% c("o", "c"))] <- 0
    y[which(y %in% c("i", ":"))] <- 1
    
    x <- which(c(-9, y) != c(y, -9))
    ifelse(1 %in% diff(x), "Singleton.Present", "Singleton.Absent")
  }
  )
  )
}

rectab$SingletonStatus <- SingletonStatusFunc(rectab$data)

ggplot(subset(rectab, SingletonStatus == "Singleton.Absent"),aes(RecombCount)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~Chr, scales = "free")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Deal with double recombinants                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab.sex <- subset(rectab, Chr == 27)
rectab <- subset(rectab, Chr != 27)

head(rectab[,-which(names(rectab) %in% c("data"))])
rectab$UniqueID <- paste(rectab$Chr, rectab$Family, rectab$RRID, rectab$parent, sep = "_")

#~~ create a smaller table to determine the recombination break points

distab <- subset(rectab, select = c(Chr, Family, RRID, parent, RecombCount, UniqueID, data, Offspring.SEX))

#~~ Now create a function to pull out the segment lengths

switchPosFunc <- function(string, UniqueID){
  
  #~~ recode phased alleles as 0 and 1 and remove uninformative loci
  
  y <- unlist(strsplit(string, split = ""))
  y[which(y %in% c("o", "c"))] <- 0
  y[which(y %in% c("i", ":"))] <- 1
  
  y <- data.frame(Phase = y)
  y$Order <- 1:nrow(y)
  
  y <- y[-which(y$Phase == "-"),]
  y$Phase <- as.numeric(as.character(y$Phase))
  
  #~~ Determine the first and last positions of each segment between crossover points
  
  y <- rbind(y, c(2, -999))
  
  y$Temp  <- c(2, y$Phase[1:(nrow(y)-1)])
  
  sortvec <- sort(c(which(y$Phase != y$Temp)))
  
  sortvec <- sort(c(sortvec, sortvec[2:length(sortvec)]-1))
  
  y2 <- y[sortvec[-length(sortvec)],-ncol(y)]
  y2$Inf.Order <- sortvec[-length(sortvec)]
  
  #~~ Create a table of the segments
  
  if(nrow(y2) > 0){
    x <- data.frame(Phase = y2$Phase[seq(1, nrow(y2), 2)],
                    StartPos = y2$Order[seq(1, nrow(y2), 2)],
                    StopPos = y2$Order[seq(2, nrow(y2), 2)],
                    StartInf = y2$Inf.Order[seq(1, nrow(y2), 2)],
                    StopInf = y2$Inf.Order[seq(2, nrow(y2), 2)])
    x$StartSpan <- c(1, x$StopPos[-nrow(x)])
    x$StopSpan  <- c(x$StartPos[-1], x$StopPos[nrow(x)])
    
    x$Segment <- 1:nrow(x)
    x$Segment.Count <- nrow(x)
    x$Type <- "Mid"
    x$Type[1] <- "First"
    x$Type[nrow(x)] <- "Last"
    if(nrow(x) == 1) x$Type <- "Only"
    x$UniqueID <- UniqueID
    
  }
  
  if(nrow(y2) == 0){
    x <- data.frame(Phase = NA,
                    StartPos = NA,
                    StopPos = NA,
                    StartInf = NA,
                    StopInf = NA,
                    StartSpan = NA,
                    StopSpan = NA,
                    Segment = NA,
                    Segment.Count = NA,
                    Type = NA,
                    UniqueID = UniqueID)
  }
  
  x
  
}  

#~~ Apply function using lapply and create a table with a unique line for each segment

system.time(test <- mapply(switchPosFunc, distab$data, distab$UniqueID, SIMPLIFY = F))

switchtab <- data.frame(data.table::rbindlist(test))

# system.time(switchtab <- do.call("rbind", test))
row.names(switchtab) <- 1:nrow(switchtab)

#~~ Merge with distab and add map information

switchtab <- join(switchtab, distab)
switchtab <- subset(switchtab, select = -data)

switchmap <- subset(lmap, select = c(Chr, Order, GenomicPosition, cMPosition))

names(switchmap) <- c("Chr", "StopSpan", "StopSpan.GenomicPosition", "StopSpan.cMPosition")
switchtab <- join(switchtab, switchmap, by = c("Chr", "StopSpan"))

names(switchmap) <- c("Chr", "StartSpan", "StartSpan.GenomicPosition", "StartSpan.cMPosition")
switchtab <- join(switchtab, switchmap, by = c("Chr", "StartSpan"))

names(switchmap) <- c("Chr", "StopPos", "Stop.GenomicPosition", "Stop.cMPosition")
switchtab <- join(switchtab, switchmap, by = c("Chr", "StopPos"))

names(switchmap) <- c("Chr", "StartPos", "Start.GenomicPosition", "Start.cMPosition")
switchtab <- join(switchtab, switchmap, by = c("Chr", "StartPos"))

#~~ Sort by Unique ID and StartPos

switchtab <- switchtab[with(switchtab, order(UniqueID, StartPos)),]
head(switchtab)

#~~ Calculate various measures related to each segment

# The genomic distance from the first to last SNP in the segment
switchtab$Segment.Length       <- switchtab$Stop.GenomicPosition - switchtab$Start.GenomicPosition
# The cM distance from the first to last SNP in the segment
switchtab$Segment.cMLength     <- switchtab$Stop.cMPosition - switchtab$Start.cMPosition
# The number of informative markers in the segment
switchtab$Segment.MarkerLength <- (switchtab$StopInf - switchtab$StartInf) + 1
# The genomic distance from the last SNP of the previous segment to the first SNP of the following segment
switchtab$Segment.SpanLength   <- switchtab$StopSpan.GenomicPosition - switchtab$StartSpan.GenomicPosition

switchtab$Singleton <- ifelse(switchtab$Segment.MarkerLength == 1, "yes", "no")


# Plot the genomic distance from the first to last SNP in the segment
ggplot(switchtab, aes(log10(Segment.Length)))     + geom_histogram(binwidth = 0.1)
# Plot the cM distance from the first to last SNP in the segment
ggplot(switchtab, aes(log10(Segment.cMLength)))   + geom_histogram(binwidth = 0.1)
# Plot the number of informative markers in the segment
ggplot(switchtab, aes(Segment.MarkerLength))      + geom_histogram() + facet_wrap(~ Type)
# Plot the genomic distance from the last SNP of the previous segment to the first SNP of the following segment
ggplot(switchtab, aes(log10(Segment.SpanLength), fill = Singleton)) + 
  geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ Type) +
  scale_fill_brewer(palette = "Set1")


#~~ Examine the X chromosome in mothers

switchtab <- join(switchtab, distab)

# Plot the genomic distance from the first to last SNP in the segment
ggplot(subset(switchtab, Chr == 27 & parent == "MOTHER" & Offspring.SEX == "Female"), aes(log10(Segment.Length)))     + 
  geom_histogram(binwidth = 0.1)
# Plot the cM distance from the first to last SNP in the segment
ggplot(subset(switchtab, Chr == 27 & parent == "MOTHER" & Offspring.SEX == "Female"), aes(log10(Segment.cMLength)))   + 
  geom_histogram(binwidth = 0.1)
# Plot the number of informative markers in the segment
ggplot(subset(switchtab, Chr == 27 & parent == "MOTHER" & Offspring.SEX == "Female"), aes(Segment.MarkerLength))      + 
  geom_histogram() + facet_wrap(~ Type)
# Plot the genomic distance from the last SNP of the previous segment to the first SNP of the following segment
ggplot(subset(switchtab, Chr == 27 & parent == "MOTHER" & Offspring.SEX == "Female"), aes(log10(Segment.SpanLength), fill = Singleton)) + 
  geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ Type) +
  scale_fill_brewer(palette = "Set1")


ggplot(subset(switchtab, Type == "Mid" & Chr %in% 1:26), aes(Segment.MarkerLength, log10(Segment.SpanLength))) + 
  geom_point(alpha = 0.1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x = "Number of SNPs", y = "log10 Span Width (Mb)")
  


nrow(subset(switchtab, Singleton == "no"))
switchtab$switchtab2
ggplot(subset(switchtab, Type == "Mid" & Chr %in% 1:26), aes(Segment.SpanLength, fill = Singleton)) +
  geom_histogram(binwidth = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Span length of double crossover (Mb)", fill = "Single SNP Crossover?") +
  geom_vline(xintercept = c(9.7e6)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "top")
  
  #geom_vline(xintercept = c(5e6), linetype = "dashed")

write.table(switchtab, "results/2_switchtab.doublexovers.forpaperwithSingletons.txt",row.names = F, quote = F, sep = "\t")



nrow(subset(switchtab, Type == "Mid" & Chr %in% 1:26 & Segment.SpanLength < 9.7e6 & Segment.SpanLength > 5e6))
(297*2)/(98420+(297*2))


median(switchtab$Segment.SpanLength)

#~~ Remove singletons

removelines <- c(which(switchtab$Singleton == "yes"))

nrow(switchtab)
removesections <- switchtab[removelines,]

switchtab <- switchtab[-removelines,]
nrow(switchtab)



switchtab$Segment.cMSpanLength <- switchtab$StopSpan.cMPosition - switchtab$StartSpan.cMPosition

ggplot(subset(switchtab, Type == "Mid"), aes(Segment.cMSpanLength)) +
  geom_histogram() +
  #defopts +
  labs(x = "Length of Segment between two crossovers (Log10 bp)")

head(lmap)

ggplot(lmap, aes(r, cMdiff)) + geom_point()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Check the Mid section to identify short recombinant fragments #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

switchtab.Mid   <- subset(switchtab, Type == "Mid" & Chr %in% c(1:26))

ggplot(switchtab.Mid, aes(log10(Segment.SpanLength))) + 
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~Chr, scales = "free")

mean.span <- mean(log10(switchtab.Mid$Segment.SpanLength))
sd.span   <- sd(log10(switchtab.Mid$Segment.SpanLength))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Create new data and recombination Count after dealing with removesections #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

RecCountFunc <- function(vector){
  
  unlist(lapply(vector, function(test){
    
    y <- unlist(strsplit(test, split = ""))
    y <- y[which(y != "-")]
    y[which(y %in% c("o", "c"))] <- 0
    y[which(y %in% c("i", ":"))] <- 1
    
    y1 <- length(which(c(-9, y) != c(y, -9))) - 2
    if(y1 == -2) y1 <- NA
    
    return(y1)
  }
  )
  )
}



tail(removesections)

rectab$data.v2 <- rectab$data

for(i in 1:nrow(removesections)){
  
  if(i %in% seq(1, nrow(removesections), 100)) print(paste("Fixing Problem", i, "of", nrow(removesections)))
  
  x <- rectab$data.v2[rectab$UniqueID == removesections$UniqueID[i]]
  x <- unlist(strsplit(x, split = ""))
  x[removesections$StartPos[i]:removesections$StopPos[i]] <- "-"
  x <- paste(x, collapse = "")
  
  
  rectab$data.v2[rectab$UniqueID == removesections$UniqueID[i]] <- x
  
}

rectab$RecombCount.v2 <- RecCountFunc(rectab$data.v2)

ggplot(rectab, aes(RecombCount, RecombCount.v2)) + geom_point(alpha = 0.1)

ggplot(rectab, aes(Chr, RecombCount.v2)) + geom_boxplot()

#~~ separate rectab into autosomes and sex chromosomes

recomb.summ <- data.frame(TotalRecombCount = tapply(rectab$RecombCount.v2, rectab$RRID.Offspring.ID, sum))
recomb.summ$RRID.Offspring.ID <- row.names(recomb.summ)
recomb.summ$RRID <- unlist(lapply(recomb.summ$RRID.Offspring.ID, function(string) unlist(strsplit(string, split = "_")[[1]][1])))
recomb.summ$Offspring.ID <- unlist(lapply(recomb.summ$RRID.Offspring.ID, function(string) unlist(strsplit(string, split = "_")[[1]][2])))

ggplot(recomb.summ, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ 4. Create new tables for further analysis                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dput(names(rectab))
rectab.v2 <- rectab[,c("Offspring.ID", "parent", "Chr", "Family",
                       "RRID", "Chromosome.Length", "Chromosome.SNP.Count", 
                       "FATHER", "MOTHER", "RRID.SEX", "RRID.BYEAR", "Offspring.SEX", "Offspring.BYEAR", 
                       "RRID.CAPAGE", "RRID.Offspring.ID", "UniqueID",
                       "data.v2", "RecombCount.v2")]

#~~ Determine the number of informative loci

NoInfLociFunc <- function(string){
  temp <- which(unlist(strsplit(string, split = "")) %in% c("i", "1","0","o", ":", "c"))
  length(temp)
}  

rectab.v2$No.Inf.Loci <- unlist(lapply(rectab.v2$data.v2, NoInfLociFunc))


#~~ Determine the Informative Start and Stop Positions

InfLengthFunc <- function(CMPstring, position = "First"){
  test <- gsub(" ", "", CMPstring)
  temp <- which(unlist(strsplit(test, split = "")) %in% c("i", "1","0","o", ":", "c"))
  if(position == "First") return(temp[1])
  if(position == "Last"){
    if(length(temp) == 0) return(NA)
    if(length(temp) > 0) return(temp[length(temp)])
  }
}  

rectab.v2$First.Inf.Order <- unlist(lapply(rectab.v2$data.v2, InfLengthFunc))
rectab.v2$Last.Inf.Order <-  unlist(lapply(rectab.v2$data.v2, function(x) InfLengthFunc(x, position = "Last")))

#~~ Add position information
head(lmap)

lmap$First.Inf.Order <- lmap$Order
lmap$Last.Inf.Order  <- lmap$Order

maptabstart <- subset(lmap, select = c(Chr, GenomicPosition, First.Inf.Order))
maptabstop  <- subset(lmap, select = c(Chr, GenomicPosition, Last.Inf.Order ))

names(maptabstart)[2] <- "First.Inf.Pos"
names(maptabstop)[2] <- "Last.Inf.Pos"

rectab.v2 <- merge(rectab.v2, maptabstart, all.x = T)
rectab.v2 <- merge(rectab.v2, maptabstop , all.x = T)

#~~ Determine the proportion of the chromosome captured between the first and last markers
rectab.v2$Prop.Inf.Chromosome <- (rectab.v2$Last.Inf.Pos - rectab.v2$First.Inf.Pos)/rectab.v2$Chromosome.Length

rectab.v2$Prop.Inf.SNPs <- rectab.v2$No.Inf.Loci/rectab.v2$Chromosome.SNP.Count

#~~ determine grandparental phase information

#~~ create temporary columns in the results table
rectab.v2$temp <- gsub(" ", "", rectab.v2$data)
rectab.v2$temp <- gsub("-", "", rectab.v2$temp)

rectab.v2$tempGFather <- gsub("o", "", rectab.v2$temp)
rectab.v2$tempGFather <- gsub("0", "", rectab.v2$tempGFather)
rectab.v2$tempGFather <- gsub(":", "", rectab.v2$tempGFather)
rectab.v2$tempGFather <- gsub("c", "", rectab.v2$tempGFather)

rectab.v2$GrandPat.Count <- nchar(rectab.v2$tempGFather)

rectab.v2$tempGMother <- gsub("i", "", rectab.v2$temp)
rectab.v2$tempGMother <- gsub("1", "", rectab.v2$tempGMother)
rectab.v2$tempGMother <- gsub(":", "", rectab.v2$tempGMother)
rectab.v2$tempGMother <- gsub("c", "", rectab.v2$tempGMother)

rectab.v2$GrandMat.Count <- nchar(rectab.v2$tempGMother)

rectab.v2 <- subset(rectab.v2, select = -c(temp, tempGFather, tempGMother))

#~~ determine which proportion of informative SNPs came from which grandparent
rectab.v2$Prop.GrandPat.SNPs <- rectab.v2$GrandPat.Count/rectab.v2$No.Inf.Loci
rectab.v2$Prop.GrandMat.SNPs <- rectab.v2$GrandMat.Count/rectab.v2$No.Inf.Loci


#~~ informative chromosome length

rectab.v2$Inf.Chr.Length <- rectab.v2$Last.Inf.Pos - rectab.v2$First.Inf.Pos

rectab.v2$RecombRate <- rectab.v2$RecombCount.v2/rectab.v2$Inf.Chr.Length

#~~ determine non-recombinant chromosomes


uninfchrcount <- function(vector) {
  x <- data.frame(table(vector))
  ifelse(x[1,1] == 0, x[1,2], 0)
}

#~~ Begin the Summary Table

rectab.v2$UniqueID2 <- paste(rectab.v2$Family, rectab.v2$RRID, sep = "_RRID")

rectab.x.v2 <- subset(rectab.v2, Chr == 27)
rectab.v2   <- subset(rectab.v2, Chr != 27)


head(rectab.v2[,-which(names(rectab.v2) == "data.v2")])

recsumm <- data.frame(
  # Sum of all recombination events including non-recombinant chromosomes
  TotalRecombCount = tapply(rectab.v2$RecombCount.v2, rectab.v2$UniqueID2, sum, na.rm = T),
  # Count of all genome-wide informative loci
  TotalInfLoci     = tapply(rectab.v2$No.Inf.Loci, rectab.v2$UniqueID2, sum),
  # Mean proportion of the chromosomes covered by informative loci
  MeanPropChr      = tapply(rectab.v2$Prop.Inf.Chromosom, rectab.v2$UniqueID2, mean),
  # Mean proportion of Loci on the chromosomes which are informative
  MeanPropInfLoci  = tapply(rectab.v2$Prop.Inf.SNPs, rectab.v2$UniqueID2, mean),
  # Number of chromosomes with NO informative loci
  UninfChrCount    = tapply(rectab.v2$No.Inf.Loci, rectab.v2$UniqueID2, uninfchrcount),    
  # How many of the chromosomes are not recombinant?
  NonRecombChrCount= tapply(rectab.v2$RecombCount.v2, rectab.v2$UniqueID2, uninfchrcount),
  # Mean number of informative loci inherited from grandfather
  MeanGrandPaternal= tapply(rectab.v2$GrandPat.Count, rectab.v2$UniqueID2, mean),
  # Mean number of informative loci inherited from grandmother
  MeanGrandMaternal= tapply(rectab.v2$GrandMat.Count, rectab.v2$UniqueID2, mean),
  # Total Informative Chr Length INCLUDING non-recombinant chromosomes
  TotalInfChrLenIncNonRecombs= tapply(rectab.v2$Inf.Chr.Length, rectab.v2$UniqueID2, function (x) sum(as.numeric(x))),
  # Mean RR including ALL chromosomes (even non-recombinants)
  MeanRRincNonRecombs           = tapply(rectab.v2$RecombRate, rectab.v2$UniqueID2, mean),
  # Variance in RR including ALL chromosomes (even non-recombinants)
  VarRRincNonRecombs            = tapply(rectab.v2$RecombRate, rectab.v2$UniqueID2, var))

#~~ add ID information

recsumm$UniqueID2 <- row.names(recsumm)
recsumm$RRID <- NA
recsumm$Family <- NA
recsumm$Offspring.ID <- NA

for(i in 1:nrow(recsumm)){
  x <- unlist(strsplit(recsumm$UniqueID2[i], split = "_RRID"))
  recsumm$RRID[i] <- x[2]
  recsumm$Family[i] <- x[1]
  recsumm$Offspring.ID[i] <- gsub("Offspring__", "", recsumm$Family[i])
  rm(x)
}


head(recsumm)

#~~ merge with sex etc.

phenoinfo <- unique(subset(rectab.v2, select = c(RRID, RRID.SEX, RRID.BYEAR, Offspring.ID, Offspring.SEX, Offspring.BYEAR, RRID.CAPAGE)))

head(phenoinfo)
head(recsumm)

recsumm$Offspring.ID <- gsub("Offspring_Dad_", "", recsumm$Offspring.ID)
recsumm$Offspring.ID <- gsub("Offspring_Mum_", "", recsumm$Offspring.ID)


recsumm <- merge(recsumm, phenoinfo, by = c("RRID", "Offspring.ID"))


ggplot(recsumm, aes(TotalRecombCount)) + 
  geom_histogram(binwidth = 1) +
  facet_wrap(~RRID.SEX)



#~~ write to file

write.table(rectab.v2,
            paste("results/2_IndivRRDoubClean_", AnalysisSuffix, ".txt", sep = ""), 
            row.names = F, sep = "\t", quote = F)

write.table(rectab.x.v2,
            paste("results/2_IndivRRDoubClean_Xchr_", AnalysisSuffix, ".txt", sep = ""), 
            row.names = F, sep = "\t", quote = F)

write.table(recsumm,
            paste("results/2_TotalSummIndivRRDoubClean_", AnalysisSuffix, ".txt", sep = ""), 
            row.names = F, sep = "\t", quote = F)













