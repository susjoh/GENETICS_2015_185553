#
#
# Characterisation of the recombination landscape
# Susan Johnston
#
#
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ load functions and libraries

library(beepr)
library(ggplot2)
library(reshape)
library(plyr)
library(GenABEL)

source("r/multiplot.R")

options(stringsAsFactors = F)

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in data

# Map data 

maptab <- read.table(paste0("results/2_Merged_map_", AnalysisSuffix, ".txt"), header = T)
tail(maptab)
maptab$Cumu <- cumsum(c(maptab$GenomicPosition[1], as.numeric(maptab$GenomicDiff)))[1:(nrow(maptab))]
maxvals <- read.table("results/2_MaxVals_g.txt", header = T)

# defopts

defopts <- theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 strip.background = element_blank())

sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")


gc.content <- read.table("C:/Users/Susan Johnston/Desktop/Recombination Rate Study/results/superceded/2.3_GC_Content_100kb_e.txt", header = T)
head(gc.content)

# Functions

binpr <- function(r, GenomicPosition, chr, binsize, windowshift = NULL){
  
  require(plyr)
  r[which(is.na(r))] <- 0

  test <- data.frame(r = r,
                     GenomicPosition = GenomicPosition,
                     Chr = chr)
  
  test <- arrange(test, Chr, GenomicPosition)
  
  test$Bin <- as.numeric(as.character(cut(test$GenomicPosition,
                                          seq(1, max(test$GenomicPosition) + binsize, binsize),
                                          labels = seq(1, max(test$GenomicPosition), binsize))))
  
  if(!is.null(windowshift)){
    test$Bin <- as.numeric(as.character(cut(test$GenomicPosition,
                                          seq(windowshift, max(test$GenomicPosition) + binsize, binsize),
                                          labels = seq(windowshift, max(test$GenomicPosition), binsize))))
    test$Bin[which(is.na(test$Bin))] <- 1
  }
  
  test$BinFirst <- c(-9, test$Bin[-length(test$Bin)])
  test$BinLast  <- c(test$Bin[-1], -9)
  
  test$Position <- ifelse(test$Bin != test$BinFirst, "First",
                          ifelse(test$Bin != test$BinLast, "Last", "Mid"))
  
  test <- subset(test, select = -c(BinFirst, BinLast))
  
  test$DistanceToBinStart <- test$GenomicPosition - test$Bin
  test$DistanceToBinEnd <- test$Bin + binsize - 1 - test$GenomicPosition
  test$DistanceToPrevSNP<- c(test$GenomicPosition[1], diff(test$GenomicPosition)) - 1
  test$DistanceToNextSNP<- c(diff(test$GenomicPosition), (test$Bin[nrow(test)] + binsize) - test$GenomicPosition[nrow(test)])
  test$Prev.r<- c(0, test$r[-length(test$r)])
  
  test$DistanceToPrevSNP[which(test$DistanceToPrevSNP < 0)] <- test$DistanceToBinStart  [which(test$DistanceToPrevSNP< 0)]
  test$DistanceToNextSNP[which(test$DistanceToNextSNP < 0)] <- test$DistanceToBinEnd  [which(test$DistanceToNextSNP < 0)]

  
  #~~ The first marker will absorb the probability of crossover from the beginning of the bin
  
  test$PrXover <- ifelse(test$Position == "First", (test$DistanceToBinStart/test$DistanceToPrevSNP * test$Prev.r) + test$r,
                         ifelse(test$Position == "Mid", test$r,
                                ifelse(test$Position == "Last", test$DistanceToBinEnd/test$DistanceToNextSNP * test$r, NA)))
  
  

  result.tab <- melt(tapply(test$PrXover, list(test$Bin, test$Chr), sum))
  names(result.tab) <- c("Bin", "Chr", "PrXover")
  result.tab <- result.tab[which(!is.na(result.tab$PrXover)),]
  
  result.tab
  
  }

model.test <- function(x, y){
  
  library(AICcmodavg); library(plyr); library(stringr)
  dat <- data.frame(x = x,
                    y = y)  
  models <- list(lm(y~x, data = dat), 
                 lm(y~I(1/x), data=dat),
                 lm(y ~ log(x), data = dat),
                 lm(y ~ I(x^2) + x, data = dat),
                 nls(y ~ I(1/x*a) + b*x, data = dat, start = list(a = 1, b = 1)), 
                 nls(y ~ (a + b*log(x)), data=dat, start = setNames(coef(lm(y ~ log(x), data=dat)), c("a", "b"))),
                 nls(y ~ I(exp(1)^(a + b * x)), data=dat, start = list(a=0,b=0)),
                 nls(y ~ I(1/x*a)+b, data=dat, start = list(a=1,b=1)))
  
  
  # have a quick look at the visual fit of these models
  print(ggplot(dat, aes(x, y)) + geom_point(alpha = 0.2) +
          stat_smooth(method = "lm", formula = as.formula(models[[1]]), size = 1, se = FALSE, colour = "black") + 
          stat_smooth(method = "lm", formula = as.formula(models[[2]]), size = 1, se = FALSE, colour = "blue") + 
          stat_smooth(method = "lm", formula = as.formula(models[[3]]), size = 1, se = FALSE, colour = "yellow") + 
          stat_smooth(method = "lm", formula = as.formula(models[[4]]), size = 1, se = FALSE, colour = "green") + 
          stat_smooth(method = "nls", formula = as.formula(models[[5]]), data=dat, start = list(a=0,b=0), size = 1, se = FALSE, colour = "red") + 
          stat_smooth(method = "nls", formula = as.formula(models[[6]]), data=dat, start = setNames(coef(lm(y ~ log(x), data=dat)), c("a", "b")), size = 1, se = FALSE, colour = "pink") +
          stat_smooth(method = "nls", formula = as.formula(models[[7]]), data=dat, start = list(a=0,b=0), size = 1, se = FALSE, colour = "violet") +
          stat_smooth(method = "nls", formula = as.formula(models[[8]]), data=dat, start = list(a=0,b=0), size = 1, se = FALSE, colour = "orange")
  )
  
  cbind(ldply(models, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) }),
        Colour = c("black", "blue", "yellow", "green", "red", "pink", "violet", "orange"))
  #ldply(models, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Broad Scale Variation                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(maptab, aes(GenomicPosition/1e6, cMPosition)) + 
  geom_point(alpha = 0.1) +
  facet_wrap(~Chr, scale = "free") +
  labs(x = "Genomic Position (Mb)", y = "Linkage Map Position (cM)") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        #strip.background = element_blank(),
        legend.position = "top") 


#~~ Create a graph of sex-specific linkage map information for plotting

maptab.sex <- melt(maptab[,c("Chr", "GenomicPosition", "cMPosition.Female", "cMPosition.Male")], id.vars = c("Chr", "GenomicPosition"))
head(maptab.sex)
maptab.sex$variable <- gsub("cMPosition.", "", maptab.sex$variable)

head(maptab.sex)
maptab.sex$Chr2 <- maptab.sex$Chr
maptab.sex$Chr2 <- gsub("27", "X", maptab.sex$Chr2)

#~~ Create a plot of the per chromosome linkage maps

#pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\PerChrLinkageMaps.pdf", width = 10, height = 8) 
ggplot(maptab.sex, aes(GenomicPosition/1e6, value, col = variable)) + 
  geom_line() +
  scale_colour_manual(values = c("black", "darkgrey")) +
  facet_wrap(~Chr, scale = "free") +
  labs(x = "Genomic Position (Mb)", y = "Linkage Map Position (cM)", col = "Sex") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 8, vjust = 0),
        axis.text.y  = element_text (size = 8, hjust = 1.3),
        #strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.background = element_blank(),
        legend.position = "top") 
#dev.off()



ggplot(maptab.sex, aes(GenomicPosition/1e6, value, col = variable)) + 
  geom_point() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Chr, scale = "free") +
  labs(x = "Genomic Position (Mb)", y = "Linkage Map Position (cM)", col = "Sex") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, vjust = 0),
        axis.text.y  = element_text (size = 12, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top") 



#~~ Get map statistics and write to LaTeX table

sum(as.numeric(maxvals$maxGenome))/nrow(maptab)
maptab$GenomicDiff[which(maptab$GenomicDiff == 0)] <- NA

median(maptab$GenomicDiff, na.rm = T)
maptab$Bin <- cut(maptab$GenomicDiff, breaks = seq(-500, max(maptab$GenomicDiff, na.rm = T), 1000))
head(arrange(data.frame(table(maptab$Bin)), -Freq))

newmaxvals <- maxvals[,c(5, 6, 1:4)]
for(i in 3:5) newmaxvals[,i]<- format(round(newmaxvals[,i], 1), nsmall = 1)
# write.table(newmaxvals, "latex/chromosomeinfo.txt", sep = " & ", quote = F, row.names = F, eol = "\\\\\n")


#~~ Create plots summarising genomc v linkage, recombination rate, and heterochiasmy

summary(lm(maxcM ~ maxGenome, data = maxvals))

maxvals$Chr2 <- maxvals$Chr
maxvals$Chr2 <- gsub("27", "X", maxvals$Chr2)

#pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\ChromosomeInfo.pdf", width = 15, height = 5) 
multiplot(
  ggplot(subset(maxvals, Chr != 27), aes(maxGenome/1e6, maxcM)) + 
    stat_smooth(method = "lm") +
    geom_text(aes(label = Chr2)) +
    # geom_text(x = maxvals$maxGenome[27]/1e6, y = maxvals$maxcM[27], label = "X") +
    labs(x = "Chromosome Length (Mb)", y = "Chromosome Length (cM)") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          plot.title = element_text(size = 18, hjust = 0),
          strip.background = element_blank()) +
    ggtitle("A") + 
    scale_x_continuous(breaks = seq(0, 300, 50)) ,
  
  ggplot(subset(maxvals, Chr != 27), aes(maxGenome/1e6, maxcM/maxGenome*1e6)) + 
    stat_smooth(method = "lm", formula = y ~ I(1/x)) +
    geom_text(aes(label = Chr2)) +
    labs(x = "Chromosome Length (Mb)", y = "Recombination Rate (cM/Mb)") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          plot.title = element_text(size = 18, hjust = 0),
          strip.background = element_blank()) +
    ggtitle("B") + 
    scale_x_continuous(breaks = seq(0, 300, 50)),
  
  ggplot(subset(maxvals, Chr != 27), aes(maxcM.male, maxcM.female)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = Chr)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          plot.title = element_text(size = 18, hjust = 0),
          strip.background = element_blank(),
          legend.position = "top") +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  ggtitle("C") + 
  labs(x = "Male Linkage Map Length (cM)", y = "Female Linkage Map Length (cM)"),

  cols = 3)
#dev.off()

#~~ Statistics!

summary(lm(maxcM ~ maxGenome, data = subset(maxvals, Chr != 27)))
summary(lm(I(maxGenome/maxcM) ~ I(1/maxcM), data = subset(maxvals, Chr != 27)))
summary(lm(maxcM.female ~ maxcM.male, data = subset(maxvals, Chr != 27)))


#~~ Sex specific

#~~ Plot the maxvals by sex

sexvals <- melt(maxvals, id.vars = c("maxGenome", "Chr", "nsnps", "Chr2"))
sexvals$variable <- gsub("maxcM.", "", sexvals$variable)
sexvals$variable <- gsub("maxcM", "Both sexes", sexvals$variable)
sexvals$value <- as.numeric(sexvals$value)

#pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\ChromosomeInfoSexSpecific.pdf", width = 12, height = 6) 
multiplot(
  ggplot(subset(sexvals, variable != "Both sexes"), aes(maxGenome/1e6, value, col = variable)) + 
    stat_smooth(method = "lm") +
    geom_text(aes(label = Chr)) +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Chromosome Length (Mb)", y = "Chromosome Length (cM)", col = "") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          plot.title = element_text(size = 18, hjust = 0, vjust = -4),
          strip.background = element_blank(),
          legend.position = "top") +
    ggtitle("A") + 
    scale_x_continuous(breaks = seq(0, 300, 50))
  ,
  
  ggplot(subset(sexvals, variable != "Both sexes"), aes(maxGenome/1e6, value/maxGenome*1e6, col = variable)) + 
    stat_smooth(method = "lm", formula = y ~ I(1/x)) +
    geom_text(aes(label = Chr)) +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Chromosome Length (Mb)", y = "Recombination Rate (cM/Mb)", col = "") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          plot.title = element_text(size = 18, hjust = 0, vjust = -4),
          strip.background = element_blank(),
          legend.position = "top") +
    ggtitle("B") + 
    scale_x_continuous(breaks = seq(0, 300, 50))
  ,cols = 2)
#dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Regional Variation in recombination          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(maptab)
head(maxvals)

maptab <- join(maptab, maxvals[,c("Chr", "maxGenome")])

#~~ Get MAF info

alleletab <- summary.snp.data(gtdata(sheepabel))
head(alleletab)
alleletab$SNP.Name <- row.names(alleletab)

#~~ Get informative meiosis information and join to alleletab

meiotab <- NULL
for(i in 1:27){
  x <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".loc"), skip = 5)
  x$Chr <- i
  meiotab <- rbind(meiotab, x)
}

meiotab <- meiotab[,-1]
names(meiotab) <- c("SNP.Name", "no.inf.mei", "no.inf.mei.PK", "tot_f", "tot_m",  "pk_f",  "pk_m", "Chromosome")

alleletab <- join(alleletab, meiotab)
head(alleletab)
  
#~~ Calculate crossover probability per bin

binsize <- 1000000

bothsex.pr <- binpr(maptab$r, maptab$GenomicPosition, maptab$Chr, binsize)
male.pr <- binpr(maptab$Male.r, maptab$GenomicPosition, maptab$Chr, binsize)
female.pr <- binpr(maptab$Female.r, maptab$GenomicPosition, maptab$Chr, binsize)

names(bothsex.pr)[3] <- "Sex.averaged.Pr"
names(male.pr)   [3] <- "Male.Pr"
names(female.pr) [3] <- "Female.Pr"

probtab <- join(bothsex.pr, male.pr)
probtab <- join(probtab,    female.pr)

#~~ Remove sex chromosomes for just now

probtab <- droplevels(subset(probtab, Chr != 27))
tapply(maptab$r, maptab$Chr, sum, na.rm = T)

#~~ Add other information on chromosome to data.frame 

probtab$BinID <- ((probtab$Bin - 1)/binsize) + 1
maxbin <- data.frame(MaxBin = tapply(probtab$BinID, probtab$Chr, max))
maxbin$Chr <- row.names(maxbin)
probtab <- join(probtab, maxbin)
probtab$DistanceToEnd <- (probtab$MaxBin - probtab$BinID) + 1
probtab$DistanceToTelo <- apply(probtab[,c("BinID", "DistanceToEnd")], 1, min)

probtab$MaleToFemale    <- probtab$Male.Pr - probtab$Female.Pr


#~~ Calculate other statistics per bin

alleletab <- droplevels(subset(alleletab, Chromosome %in% c(1:26)))
alleletab$Bin <- as.numeric(as.character(cut(alleletab$Position,
                                        seq(1, max(alleletab$Position) + binsize, binsize),
                                        labels = seq(1, max(alleletab$Position), binsize))))
alleletab$BinID <- ((alleletab$Bin - 1)/binsize) + 1

head(alleletab)

genome.stats <- melt(tapply(alleletab$Position, list(alleletab$BinID, alleletab$Chromosome), length))
names(genome.stats) <- c("BinID", "Chr", "SNP.Count")

test <- melt(tapply(alleletab$Q.2, list(alleletab$BinID, alleletab$Chromosome), mean))
names(test) <- c("BinID", "Chr", "Mean.MAF")
genome.stats <- join(genome.stats, test)

test <- melt(tapply(alleletab$no.inf.mei.PK, list(alleletab$BinID, alleletab$Chromosome), mean))
names(test) <- c("BinID", "Chr", "Mean.Inf.Mei")
genome.stats <- join(genome.stats, test)


genome.stats <- subset(genome.stats, !is.na(SNP.Count))

#~~ Check they match bin categories across tables

probtab$BinID <- as.integer(probtab$BinID)
identical(arrange(genome.stats[,c("BinID", "Chr")], BinID, Chr),arrange(probtab[,c("BinID", "Chr")], BinID, Chr))


probtab <- join(probtab, genome.stats)


#~~ gc.content

head(gc.content)
gc.content$Category <- cut(gc.content$SeqStart,
                           breaks = seq(0.5, max(gc.content$SeqStart) + binsize, binsize),
                           labels = seq(1, max(gc.content$SeqStart), binsize))

gc.content$BinID <- ((as.numeric(as.character(gc.content$Category)) - 1)/binsize) + 1

tail(gc.content)

new.gc <- melt(tapply(gc.content$Count.A, list(gc.content$BinID, gc.content$Chr), sum, na.rm = T))
names(new.gc)[3] <- "Count.A"
for(i in c("Count.C", "Count.G", "Count.T", "Count.ACGT")){
  x <- melt(tapply(gc.content[,i], list(gc.content$BinID, gc.content$Chr), sum, na.rm = T))
  names(x)[3] <- i
  new.gc <- join(new.gc, x)
}
names(new.gc)[1:2] <- c("BinID", "Chr")
new.gc <- subset(new.gc, !is.na(Count.ACGT))

new.gc$PC.GC <- rowSums(new.gc[,c("Count.C", "Count.G")])/new.gc$Count.ACGT
head(new.gc)

probtab <- join(probtab, new.gc)

#~~ Get rid of final bins if they span less than half of the binsize
test <- join(probtab[which(probtab$DistanceToEnd == 1),], maxvals)
test$LastBinSize <- test$maxGenome - test$Bin
plot(test$SNP.Count, test$LastBinSize)
test <- subset(test, LastBinSize < binsize/2)
head(test)




probtab <- probtab[-which(paste(probtab$Bin, probtab$Chr) %in% paste(test$Bin, test$Chr)),]

rm(test, female.pr, male.pr, maxbin, x, sexvals, bothsex.pr)

#~~ Divide the distance to the telomere by the length of the chromosome
probtab$AdjustedDistance <- probtab$DistanceToTelo/probtab$MaxBin

#~~ Classify telomeres based on metacentric and acrocentric. According to Jill Maddox, acrocentric
#   chromosomes should have their centromere close the beginning of the chromosome.

probtab$TelomereClass <- "noncentro"

probtab$TelomereClass[which(probtab$Chr %in% 4:26 & probtab$BinID == probtab$DistanceToTelo)] <- "centro"


#~~ Some plots for information

multiplot(
ggplot(probtab, aes(SNP.Count)) + geom_histogram(binwidth = 1, col = "grey")
,ggplot(probtab, aes(Sex.averaged.Pr)) + geom_histogram(binwidth = 0.001, col = "grey")
,ggplot(probtab, aes(Male.Pr)) + geom_histogram(binwidth = 0.001, col = "grey")
,ggplot(probtab, aes(Female.Pr)) + geom_histogram(binwidth = 0.001, col = "grey")
,ggplot(probtab, aes(MaleToFemale)) + geom_histogram(binwidth = 0.001, col = "grey") + geom_vline(xintercept = mean(probtab$MaleToFemale), colour = "red")
,ggplot(probtab, aes(Mean.MAF , Mean.Inf.Mei)) + geom_point(alpha = 0.3)
,ggplot(probtab, aes(SNP.Count, Sex.averaged.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(SNP.Count, Male.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(SNP.Count, Female.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.Inf.Mei, Sex.averaged.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.Inf.Mei, Male.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.Inf.Mei, Female.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.MAF, Sex.averaged.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.MAF, Male.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
,ggplot(probtab, aes(Mean.MAF, Female.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
, cols = 5)

ggplot(probtab, aes(PC.GC, Sex.averaged.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
ggplot(probtab, aes(PC.GC, Male.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
ggplot(probtab, aes(PC.GC, Female.Pr)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")
ggplot(probtab, aes(PC.GC, MaleToFemale)) + geom_point(alpha = 0.3) + stat_smooth(method = "lm")


ggplot(probtab, aes(Male.Pr, Female.Pr)) + geom_point(alpha = 0.2)

#~~ How does the recombination probability look per bin across the genome?

ggplot(probtab, aes(BinID, Sex.averaged.Pr)) + geom_point(alpha = 0.2) + stat_smooth(method = "loess") + facet_wrap(~Chr, scales = "free_x")
ggplot(probtab, aes(BinID, Male.Pr))         + geom_point(alpha = 0.2) + stat_smooth(method = "loess") + facet_wrap(~Chr, scales = "free_x")
ggplot(probtab, aes(BinID, Female.Pr))       + geom_point(alpha = 0.2) + stat_smooth(method = "loess") + facet_wrap(~Chr, scales = "free_x")
ggplot(probtab, aes(BinID, MaleToFemale))    + geom_point() + stat_smooth(method = "loess") + facet_wrap(~Chr, scales = "free_x")

ggplot(probtab, aes(DistanceToTelo, MaleToFemale)) + geom_point() + stat_smooth() +facet_wrap(~Chr, scales = "free_x")

ggplot(probtab, aes(BinID, PC.GC))    + geom_point() + stat_smooth(method = "loess") + facet_wrap(~Chr, scales = "free_x")

#~~ What is the interaction between the sexes?

probtab2 <- droplevels(subset(probtab, select = c(Male.Pr, Female.Pr, DistanceToTelo, AdjustedDistance, MaxBin, SNP.Count, BinID, Chr, PC.GC, TelomereClass)))
probtab2 <- melt(probtab2, id.vars = c("DistanceToTelo", "MaxBin", "SNP.Count", "BinID", "Chr", "AdjustedDistance", "PC.GC", "TelomereClass"))
head(probtab2)
names(probtab2)[which(names(probtab2) %in% c("variable", "value"))] <- c("Sex", "Probability")
probtab2$Sex <- gsub(".Pr", "", probtab2$Sex)

head(probtab2)

model.test(subset(probtab2, Sex == "Male")$DistanceToTelo, subset(probtab2, Sex == "Male")$Probability)
model.test(subset(probtab2, Sex == "Female")$DistanceToTelo, subset(probtab2, Sex == "Female")$Probability)
model.test(probtab$DistanceToTelo, probtab$MaleToFemale)

# Quadratic, multiplicative inverse and log are the best models, respectively...

ggplot(probtab, aes(DistanceToTelo, MaleToFemale)) + 
  geom_point(alpha = 0.2) + 
  #stat_smooth(method = "lm", formula = y ~ I(x^2) + x, colour = "red") + 
  #stat_smooth(method = "lm", formula =  y ~ I(1/x), colour = "blue") + 
  stat_smooth(method = "lm", formula = y ~ log10(x)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Distance to Nearest Telomere (MB)", y = "Male:Female Crossover Probability", col = "") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        strip.background = element_blank(),
        legend.position = "top")


#~~ Run some models with sex interaction

probtab4 <- subset(probtab2, DistanceToTelo < 61)

summary(glm(Probability ~ Sex * DistanceToTelo           + SNP.Count + MaxBin + PC.GC , data = probtab4, family = gaussian()))$aic
summary(glm(Probability ~ Sex * I(1/DistanceToTelo)      + SNP.Count + MaxBin + PC.GC , data = probtab4, family = gaussian()))$aic
summary(glm(Probability ~ Sex * I(log10(DistanceToTelo)) + SNP.Count + MaxBin + PC.GC , data = probtab4, family = gaussian()))$aic
summary(glm(Probability ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo 
            + SNP.Count + MaxBin + PC.GC , data = probtab4, family = gaussian()))$aic
summary(glm(Probability ~ Sex * I(DistanceToTelo^3) + Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo 
            + SNP.Count + MaxBin + PC.GC , data = probtab4, family = gaussian()))$aic

summary(glm(Probability ~ Sex * I(DistanceToTelo^3) + Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo 
            + SNP.Count + PC.GC , data = probtab4, family = gaussian()))

fit1 <- glm(Probability ~ SNP.Count + PC.GC, data = probtab4, family = gaussian())
probtab4$residuals <- residuals(fit1)

# model.test(subset(probtab4, Sex == "Male")$DistanceToTelo,
#            subset(probtab4, Sex == "Male")$residuals)
# model.test(subset(probtab4, Sex == "Female")$DistanceToTelo,
#            subset(probtab4, Sex == "Female")$residuals)



probtab3 <- subset(probtab, DistanceToTelo < 61)

summary(glm(MaleToFemale ~ DistanceToTelo + MaxBin + SNP.Count + PC.GC, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(1/DistanceToTelo) + MaxBin + SNP.Count + PC.GC, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count + PC.GC, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count + PC.GC, data = probtab3, family = gaussian()))$aic

fit2 <- glm(MaleToFemale ~ MaxBin + SNP.Count + PC.GC, data = probtab3, family = gaussian())
probtab3$residuals <- residuals(fit2)


summary(glm(residuals ~ DistanceToTelo, data = probtab3, family = gaussian()))$aic
summary(glm(residuals ~ I(1/DistanceToTelo), data = probtab3, family = gaussian()))$aic
summary(glm(residuals ~ I(log10(DistanceToTelo)), data = probtab3, family = gaussian()))$aic
summary(glm(residuals ~ I(DistanceToTelo^2) + DistanceToTelo, data = probtab3, family = gaussian()))$aic
x <- summary(glm(residuals ~ I(DistanceToTelo^3) + I(DistanceToTelo^2) + DistanceToTelo, data = probtab3, family = gaussian()))$coefficients
x <- data.frame(x)
x$Estimate


write.table(probtab, paste0("results/2_PopWide_CrossoverPr_bin", binsize/1000, "K_Full_", AnalysisSuffix, ".txt"),
            quote = F, sep = "\t", row.names = F)


model.test(probtab3$DistanceToTelo, probtab3$residuals)

#pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\CrossoverProbability.pdf", width = 9, height = 4) 
multiplot(
ggplot(probtab4, aes(DistanceToTelo, residuals, col = Sex)) +
  geom_point(alpha = 0.2) +
  #stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = "Distance to Nearest Telomere (MB)", y = "Crossover Probability", col = "") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.92))+
  ggtitle("A")
,

ggplot(probtab3, aes(DistanceToTelo, residuals)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = "Distance to Nearest Telomere (MB)", y = "Male:Female Crossover Probability", col = "") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0),
        strip.background = element_blank(),
        legend.position = c(0.8, 0.9)) +
  ggtitle("B")
, cols = 2)
#dev.off()



model2a <- summary(lm(Probability ~ Sex * I(DistanceToTelo^3) + Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count + PC.GC,
            data = probtab4,
            family = gaussian()))

model2b <- summary(lm(MaleToFemale ~ I(DistanceToTelo^3) + I(DistanceToTelo^2) + DistanceToTelo + SNP.Count + PC.GC,
            data = probtab3, family = gaussian()))

model2c <- lm(residuals ~ I(DistanceToTelo^3) + I(DistanceToTelo^2) + DistanceToTelo,
                       data = probtab3, family = gaussian())
model2a


tab1 <- data.frame(model2a$coefficients)
tab2 <- data.frame(model2b$coefficients)
tab3 <- data.frame(model2c$coefficients)


tab1
tab2

z <- matrix(c(tab3[1,1],tab3[4,1],tab3[3,1],tab3[2,1]), ncol=1) 

solve(polynomial(z))
polyroot(z)



plot(probtab3$DistanceToTelo, probtab3$residuals, pch = ".")
abline(a = 0, b = 0)
lines(probtab3$DistanceToTelo, predict(model2c), col="blue") 
abline(v = 13.1225)

install.packages("polynom")
library(polynom) 
p <- polynomial(tab3$Estimate[c(2, 3, 4, 1)]) 
solve(p) 
polyroot(tab3$Estimate[c(2, 3, 4, 1)])

# for(i in 1:4) tab1[,i]<- format(round(tab1[,i], digits = 4))

format(tab1[,1], width=7)

?prettyNum

tab1 <- cbind(Effect = row.names(tab1), tab1)
tab2 <- cbind(Effect = row.names(tab2), tab2)

newtab <- rbind(tab1, tab2)
newtab <- newtab[c(),]

for(i in 2:5) newtab[,i] <- format(newtab[,i], digits = 4)

newtab$Effect <- c("Intercept (Female)", "Male", "dist_{Telo}^{3}", "dist_{Telo}^{2}", "dist_{Telo}", 
  "SNP Count",  "Male: dist_{Telo}^{3}", "Male: dist_{Telo}^{2}", "Male: dist_{Telo}",
  "Intercept", "dist_{Telo}^{3}", "dist_{Telo}^{2}", "dist_{Telo}",
  "SNP Count") 
names(newtab) <- c("Effect", "Estimate", "S.E.", "t-value", "P")

write.table(newtab,
            file = "latex/crossover probabilities.txt",
            sep = "\t&\t", quote = F, row.names = F, eol = "\\\\\n")





# multiplot(
# ggplot(probtab2, aes(DistanceToTelo, Probability, col = Sex)) +
#   geom_point(alpha = 0.2) +
#   #stat_smooth(method = "lm", formula = y ~ log(x)) +
#   stat_smooth(method = "lm", formula = y ~ I(x^2) + x) +
#   scale_colour_brewer(palette = "Set1") +
#   labs(x = "Distance to Nearest Telomere (MB)", y = "Crossover Probability", col = "") +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 16, vjust = 0),
#         axis.text.y  = element_text (size = 14, hjust = 1.3),
#         strip.text.x = element_text (size = 16, vjust = 0.7),
#         axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
#         axis.title.x = element_text (size = 16, vjust = 0.2),
#         plot.title = element_text(size = 18, hjust = 0, vjust = -4),
#         strip.background = element_blank(),
#         legend.position = "top")
# 
# ,
# ggplot(probtab2, aes(AdjustedDistance, Probability, col = Sex)) +
#   geom_point(alpha = 0.2) +
#   #stat_smooth(method = "lm", formula = y ~ log(x)) +
#   stat_smooth(method = "lm", formula = y ~ I(x^2) + x) +
#   scale_colour_brewer(palette = "Set1") +
#   labs(x = "Distance to Nearest Telomere (MB)/\nChromosomeLength(MB)", y = "Crossover Probability", col = "") +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 16, vjust = 0),
#         axis.text.y  = element_text (size = 14, hjust = 1.3),
#         strip.text.x = element_text (size = 16, vjust = 0.7),
#         axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
#         axis.title.x = element_text (size = 16, vjust = 0.2),
#         plot.title = element_text(size = 18, hjust = 0, vjust = -4),
#         strip.background = element_blank(),
#         legend.position = "top")
# 
# ,cols = 2)



summary(glm(Sex.averaged.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic

ggplot(probtab, aes(AdjustedDistance, Sex.averaged.Pr)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(Male.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic

ggplot(probtab, aes(AdjustedDistance, Male.Pr)) +
  geom_point(alpha = 0.2) +
  #stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(Female.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab, family = gaussian()))$aic

ggplot(probtab, aes(AdjustedDistance, Female.Pr)) +
  geom_point(alpha = 0.2) +
  #stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(1/x), col = "red")




ggplot(probtab4, aes(DistanceToTelo, Probability)) +
  geom_point(alpha = 0.2) +
  #stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Distance to Nearest Telomere (MB)", y = "Crossover Probability", col = "") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        strip.background = element_blank(),
        legend.position = "top") +
  facet_wrap(~Sex)

summary(glm(Probability*100 ~ Sex * DistanceToTelo + SNP.Count, data = probtab4, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(1/DistanceToTelo) + SNP.Count, data = probtab4, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(log10(DistanceToTelo)) + SNP.Count, data = probtab4, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count, data = probtab4, family = gaussian()))$aic

summary(glm(Probability*100 ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count, data = probtab4, family = gaussian()))

# library(lme4)
# 
# summary(glmer(Probability*100 ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count + (1|Chr), data = probtab2))
# summary(glmer(Probability*100 ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count + (1|Chr), data = probtab4))



ggplot(probtab3, aes(DistanceToTelo, MaleToFemale)) + 
  geom_point(alpha = 0.2) + 
  #stat_smooth(method = "lm", formula = y ~ I(x^2) + x, colour = "red") + 
  #stat_smooth(method = "lm", formula =  y ~ I(1/x), colour = "blue") + 
  stat_smooth(method = "lm", formula = y ~ log10(x)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Distance to Nearest Telomere (MB)", y = "Male:Female Crossover Probability", col = "") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        strip.background = element_blank(),
        legend.position = "top")




summary(glm(Sex.averaged.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Sex.averaged.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic

ggplot(probtab3, aes(DistanceToTelo, Sex.averaged.Pr)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(Male.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Male.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic

ggplot(probtab3, aes(DistanceToTelo, Male.Pr)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm") +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(Female.Pr ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(Female.Pr ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic

ggplot(probtab3, aes(DistanceToTelo, Female.Pr)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(MaleToFemale ~ DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(1/DistanceToTelo) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(log10(DistanceToTelo)) + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic
summary(glm(MaleToFemale ~ I(DistanceToTelo^2) + DistanceToTelo + MaxBin + SNP.Count, data = probtab3, family = gaussian()))$aic

ggplot(probtab3, aes(DistanceToTelo, MaleToFemale)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x, col = "red")


summary(glm(Probability*100 ~ Sex * DistanceToTelo + SNP.Count, data = probtab32, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(1/DistanceToTelo) + SNP.Count, data = probtab32, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(log10(DistanceToTelo)) + SNP.Count, data = probtab32, family = gaussian()))$aic
summary(glm(Probability*100 ~ Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count, data = probtab32, family = gaussian()))$aic

ggplot(probtab32, aes(AdjustedDistance, Probability, col = Sex)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  stat_smooth(method = "lm", formula = y ~ I(x^2) + x) +
  scale_colour_brewer(palette = "Set1")







