# Create the plots of haplotypes, LD etc.
# Author: Susan Johnston

setwd("mach")

library(GenABEL)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

#~~ make dat file
source("r/recoderFunc.R")
source("../3.0.2_Beagle_Imputation_Functions.R")

# system("cmd", input = "plink --file chr6segmentmach --from s15515.1 --to s74824.1 --no-pheno --sheep --recode --out chr6segmentmachsmall")
# system("cmd", input = "plink --file chr6HDsegmentmach --from s15515.1 --no-pheno --to oar3_OAR6_117013262 --sheep --recode --out chr6HDsegmentmachsmall")

imputed <- read.table("imputed.mlgeno")[,-2]
imputed.pr <- read.table("imputed.mldose")[,-2]
imputed.snplist <- read.table("imputed.mlinfo", header = T, stringsAsFactors = F)
imputed.haplos <- read.table("snpHDset.haplos")

head(imputed)
head(imputed.pr)
head(imputed.snplist)
head(imputed.haplos)


hist(imputed.snplist$Quality)
hist(imputed.snplist$Rsq)

names(imputed) <- c("Link", as.character(imputed.snplist$SNP))
names(imputed.pr) <- c("Link", as.character(imputed.snplist$SNP))


#~~ LD plot

sheepHD <- load.gwaa.data(phe = "../../Soay Sheep HD SNP Chip/data_from_WTCRF/E10905_SheepHD_Plates1-2_310114 QC1/20140214_SheepHD_QC1_GenABELpheno.txt",
                          gen = "../../Soay Sheep HD SNP Chip/data_from_WTCRF/E10905_SheepHD_Plates1-2_310114 QC1/20140214_SheepHD_QC1_GenABEL.txt")
impabel <- load.gwaa.data(genofile = "imputedPLINK.genabel", phenofile = "imputedPLINK.phenofile")

subHD <- sheepHD[,which(chromosome(sheepHD) == 6 & map(sheepHD) > 115716891)]
nsnps(subHD)
nsnps(impabel)
subHD <- subHD[,snpnames(impabel)]


genoHD    <- data.frame(as.character.gwaa.data(subHD))
genoHDimp <- data.frame(as.character.gwaa.data(impabel))



for(i in 1:ncol(genoHD)) genoHD[,i] <- as.numeric(as.factor(genoHD[,i] ))
for(i in 1:ncol(genoHDimp)) genoHDimp[,i] <- as.numeric(as.factor(genoHDimp[,i] ))


ldtab <- data.frame(Allele1 = rep(1:ncol(genoHD), times = ncol(genoHD)),
                    Allele2 = rep(1:ncol(genoHD), each  = ncol(genoHD)),
                    HDcor = NA,
                    HDimpcor = NA)
ldtab <- ldtab[which(ldtab$Allele1 > ldtab$Allele2),]

for(i in 1:nrow(ldtab)){
  if(i %in% seq(1, nrow(ldtab), 1000)) print(i)
  
  x <- genoHD[,c(ldtab$Allele1[i],ldtab$Allele2[i])]
  x <- na.omit(x)
  
  ldtab$HDcor[i] <- cor(x[,1], x[,2])
  rm(x)
  
  x <- genoHDimp[,c(ldtab$Allele1[i],ldtab$Allele2[i])]
  x <- na.omit(x)
  
  ldtab$HDimpcor[i] <- cor(x[,1], x[,2])
  rm(x)
}

ldtab$HDcor <- ifelse(ldtab$HDcor < 0, ldtab$HDcor * -1, ldtab$HDcor)
ldtab$HDimpcor <- ifelse(ldtab$HDimpcor < 0, ldtab$HDimpcor * -1, ldtab$HDimpcor)
ldtab$Allele1.Name <- recodeChromosome


ggplot(ldtab, aes(Allele1, Allele2, fill = HDcor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ldtab$Allele1.Name <- recoderFunc(ldtab$Allele1, 1:ncol(genoHD), names(genoHD))
ldtab$Allele2.Name <- recoderFunc(ldtab$Allele2, 1:ncol(genoHD), names(genoHD))

head(ldtab)

ldtab[which(ldtab$Allele2.Name == "oar3_OAR6_116402578" & ldtab$Allele1.Name == "s74824.1"),]


#~~ Redo for the different genotypes

genoHD.AA <- subset(data.frame(as.character.gwaa.data(subHD)), oar3_OAR6_116402578 == "A/A")
genoHD.GA <- subset(data.frame(as.character.gwaa.data(subHD)), oar3_OAR6_116402578 == "G/A")
genoHD.GG <- subset(data.frame(as.character.gwaa.data(subHD)), oar3_OAR6_116402578 == "G/G")


for(i in 1:ncol(genoHD.AA)) genoHD.AA[,i] <- as.numeric(as.factor(genoHD.AA[,i] ))
for(i in 1:ncol(genoHD.GA)) genoHD.GA[,i] <- as.numeric(as.factor(genoHD.GA[,i] ))
for(i in 1:ncol(genoHD.GG)) genoHD.GG[,i] <- as.numeric(as.factor(genoHD.GG[,i] ))


ldtab.focus <- data.frame(Allele1 = rep(1:ncol(genoHD.AA), times = ncol(genoHD.AA)),
                          Allele2 = rep(1:ncol(genoHD.AA), each  = ncol(genoHD.AA)),
                          AA = NA,
                          GA = NA,
                          GG = NA)


ldtab.focus <- ldtab.focus[which(ldtab.focus$Allele1 > ldtab.focus$Allele2),]

for(i in 1:nrow(ldtab.focus)){
  if(i %in% seq(1, nrow(ldtab.focus), 1000)) print(i)
  
  x <- genoHD.AA[,c(ldtab.focus$Allele1[i],ldtab.focus$Allele2[i])]
  x <- na.omit(x)
  
  ldtab.focus$AA[i] <- cor(x[,1], x[,2])
  rm(x)
  
  x <- genoHD.GA[,c(ldtab.focus$Allele1[i],ldtab.focus$Allele2[i])]
  x <- na.omit(x)
  
  ldtab.focus$GA[i] <- cor(x[,1], x[,2])
  rm(x)
  
  x <- genoHD.GG[,c(ldtab.focus$Allele1[i],ldtab.focus$Allele2[i])]
  x <- na.omit(x)
  
  ldtab.focus$GG[i] <- cor(x[,1], x[,2])
  rm(x)
  
}

ldtab.focus$AA <- ifelse(ldtab.focus$AA < 0, ldtab.focus$AA * -1, ldtab.focus$AA)
ldtab.focus$GA <- ifelse(ldtab.focus$GA < 0, ldtab.focus$GA * -1, ldtab.focus$GA)
ldtab.focus$GG <- ifelse(ldtab.focus$GG < 0, ldtab.focus$GG * -1, ldtab.focus$GG)

ldtab.focus <- join(ldtab.focus, ldtab[,1:3])


ldtab.focus.melt <- melt(ldtab.focus, id.vars = c("Allele1", "Allele2"))

ggplot(ldtab.focus.melt, aes(Allele1, Allele2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~variable)


ggplot(ldtab.focus, aes(Allele1, Allele2, fill = AA)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ggplot(ldtab.focus, aes(Allele1, Allele2, fill = GA)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ggplot(ldtab.focus, aes(Allele1, Allele2, fill = GG)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")


#~~ :Look at haplotypes


head(imputed.haplos)

table(imputed.haplos$V3)
imputed.haplos$V3 <- toupper(imputed.haplos$V3)

imputed.haplos
imputed.haplos$V4 <- unlist(sapply(imputed.haplos$V3, function(x) substr(as.character(x), 1, 126)[[1]]))

V4tab <- data.frame(table(imputed.haplos$V4))
head(V4tab)
V4tab$HaploGroup <- letters[1:nrow(V4tab)]
head(imputed.haplos)
names(V4tab)[1] <- "V4"

imputed.haplos <- join(imputed.haplos, V4tab)

#~~~~~~~~


testfreq <- imputed.snplist[,1:4]
head(testfreq)
testfreq$refallele <- ifelse(testfreq$Freq1 <= 0.5, testfreq$Al2, testfreq$Al1)

#~~ most common haplotype

testfreq$CommonAllele <- NA
testfreq$CommonAllele[1:126] <-   strsplit(as.character(V4tab$V4[which(V4tab$Freq == max(V4tab$Freq))]), split = "")[[1]]



head(imputed.haplos)


plothaplo <- NULL

for(i in 1:nrow(imputed.haplos)){
  x <- data.frame(Allele = strsplit(as.character(imputed.haplos$V4[i]), split = "")[[1]],
                  SNP.Name = imputed.snplist$SNP[1:126],
                  HaploID = i,
                  HaploGroup = imputed.haplos$HaploGroup[i])
  x$Order <- 1:nrow(x)
  plothaplo <- rbind(plothaplo, x)
  rm(x)
}

plothaplo$Allele <- toupper(plothaplo$Allele)

head(plothaplo)




head(imputed.snplist)

test <- vector2Dendrogram(imputed.haplos$V4)$haplotype.counts
test.hc <- vector2Dendrogram(imputed.haplos$V4)$hc.object$order
test$Order <- factor(test$Order, levels = test.hc)
test <- arrange(test, Order)

head(V4tab)
V4tab$Haplo <- V4tab$V4
test <- join(test, V4tab[,c("Haplo", "HaploGroup")])

test$RNFallele <- NA
for(i in 1:nrow(test)) test$RNFallele[i] <- substr(test$Haplo[i], 66, 66)


plothaplo$FillColour <- NA
for(i in 1:nrow(plothaplo)){
  plothaplo$FillColour[i] <- ifelse(plothaplo$Allele[i] == testfreq$CommonAllele[which(testfreq$SNP == plothaplo$SNP.Name[i])], "A", "B")
}

plothaplo$HaploGroup <- factor(plothaplo$HaploGroup, levels = test$HaploGroup)
plothaplo <- arrange(plothaplo, HaploGroup)
plothaplo$NewOrder <- plothaplo$NewOrder <- rep(1:376, each = 126)

ggplot(plothaplo, aes(Order, NewOrder, fill = FillColour)) +
  geom_tile(colour = "black") +
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = 66) 

test$HaploID <- 1:nrow(test)
plothaplo.simple <- NULL

for(i in 1:nrow(test)){
  x <- data.frame(Allele = strsplit(as.character(test$Haplo[i]), split = "")[[1]],
                  SNP.Name = imputed.snplist$SNP[1:126],
                  HaploID = test$HaploID[i],
                  HaploGroup = test$HaploGroup[i],
                  HaploFreq = test$Freq[i],
                  RNFallele = test$RNFallele[i])
  x$Order <- 1:nrow(x)
  plothaplo.simple <- rbind(plothaplo.simple, x)
  rm(x)
}

plothaplo.simple$FillColour <- NA
for(i in 1:nrow(plothaplo.simple)){
  plothaplo.simple$FillColour[i] <- ifelse(plothaplo.simple$Allele[i] == testfreq$CommonAllele[which(testfreq$SNP == plothaplo.simple$SNP.Name[i])], "A", "B")
}

head(plothaplo.simple)

ggplot(plothaplo.simple, aes(Order, -HaploID, fill = FillColour)) +
  annotate("text", x = -3, y = c(0.5, -test$HaploID), label = c("Frequency", test$Freq)) +
  geom_vline(xintercept = 66) + 
  geom_tile(colour = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  scale_x_continuous(breaks = plothaplo.simple$Order[1:126],
                     labels = plothaplo.simple$SNP.Name[1:126]) +
  theme(axis.text.x  = element_text (angle = 270, vjust = 0, hjust = 0),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank (),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.margin = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  coord_cartesian(ylim = c(-14.5, 1)) +
  labs(x = "", y = "", fill = "Reference")

head(plothaplo.simple)




weehaplo <- filter(plothaplo, SNP.Name %in% c("s61186.1", "s30544.1", "s23033.1", "s14138.1", "s10844.1", "s74824.1", "oar3_OAR6_116402578"))
# 
# s61186.1 115982022 A G
# s30544.1 115991063 A G
# s23033.1 116324298 G A
# s14138.1 116441474 A G
# s10844.1 116468415 G A
# s74824.1 116668852 G A

weehaplo$Order <- as.numeric(recoderFunc(weehaplo$SNP.Name, c("s61186.1", "s30544.1", "s23033.1","oar3_OAR6_116402578", "s14138.1", "s10844.1", "s74824.1"), c(1:7)))

head(weehaplo)

weehaplo2 <- dplyr::select(weehaplo, Allele, SNP.Name, HaploID)
weehaplo2 <- cast(weehaplo2, HaploID ~ SNP.Name, value = "Allele")
head(weehaplo2)
weehaplo2$Haplotype <- apply(dplyr::select(weehaplo2, s61186.1, s30544.1, s23033.1, oar3_OAR6_116402578, s14138.1, s10844.1, s74824.1), 1, function(x) paste(x, collapse = ""))

weetest <- data.frame(table(weehaplo2$Haplotype))
weetest$NewHaplo <- 1:nrow(weetest)
names(weetest) <- c("Haplotype", "Frequency", "NewHaplo")

weeplothaplo <- NULL

for(i in 1:nrow(weetest)){
  x <- data.frame(Allele = strsplit(as.character(weetest$Haplotype[i]), split = "")[[1]],
                  SNP.Name = c("s61186.1", "s30544.1","s23033.1","oar3_OAR6_116402578", "s14138.1", "s10844.1", "s74824.1"),
                  NewHaplo = weetest$NewHaplo[i],
                  HaploFreq = weetest$Frequency[i])
  x$Order <- 1:nrow(x)
  weeplothaplo <- rbind(weeplothaplo, x)
  rm(x)
} 

weeplothaplo$NewHaplo2 <- weeplothaplo$NewHaplo
weeplothaplo$NewHaplo2[which(weeplothaplo$NewHaplo == 8)] <- 5
weeplothaplo$NewHaplo2[which(weeplothaplo$NewHaplo == 5)] <- 8


pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\SNP50haplotypes.pdf", width = 7, height = 8, useDingbats = F) 
ggplot(weeplothaplo, aes(Order, NewHaplo2, fill = Allele)) +
  geom_tile(colour = "black", alpha = 0.7) +
  geom_text(aes(label = Allele), size = 6) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  scale_x_continuous(breaks = weeplothaplo$Order[1:7],
                     labels = weeplothaplo$SNP.Name[1:7]) +
  scale_y_continuous(breaks = 1:8, labels = weetest$Frequency[c(1, 2, 3, 4, 8, 6, 7, 5)]) +
  theme(axis.text.x  = element_text (size = 12, angle = 270, vjust = 0, hjust = 0),
        axis.text.y  = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text (size = 16, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.margin = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  coord_cartesian(xlim = c(0.35, 8)) +
  labs(x = "", y = "Frequency", fill = "Allele")
dev.off()




plothaplo.simple.common <- subset(plothaplo.simple, HaploFreq > 9)
test.common <- subset(test, Freq > 9)
test.common$HaploID <- 1:nrow(test.common)
test.common

head(plothaplo.simple.common)
plothaplo.simple.common <- select(plothaplo.simple.common, -HaploID)
plothaplo.simple.common <- join(plothaplo.simple.common, select(test.common,s HaploGroup, HaploID))


ggplot(plothaplo.simple.common, aes(Order, -HaploID, fill = FillColour)) +
  annotate("text", x = -3, y = c(0.5, -test.common$HaploID), label = c("Frequency", test.common$Freq)) +
  geom_vline(xintercept = 22) + 
  geom_tile(colour = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  scale_x_continuous(breaks = plothaplo.simple.common$Order[1:82],
                     labels = plothaplo.simple.common$SNP.Name[1:82]) +
  scale_fill_manual (values = c("grey", "white")) +
  theme(axis.text.x  = element_text (angle = 270, vjust = 0, hjust = 0),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank (),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.margin = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  coord_cartesian(ylim = c(-8.5, 1)) +
  labs(x = "", y = "", fill = "Reference")











# 
# plothaplo$HaploGroup <- factor(plothaplo$HaploGroup, levels = rev(c("l", "n", "c", "j", "g", "f", "b", "m", "h", "i", "e", "k", "d", "a")))
# plothaplo <- arrange(plothaplo, HaploGroup)
# plothaplo$NewOrder <- rep(1:376, each =  82)


ggplot(plothaplo, aes(Order, HaploID, fill = FillColour)) +
  geom_tile(colour = "black") +
  scale_fill_brewer(palette = "Set1")

c("a","d","k","e","m","i","h","b","n","l","f","c","g", "j")


#















head(imputed.haplos)

table(imputed.haplos$V3)

newhaplo <- data.frame(table(imputed.haplos$V3))
newhaplo$HaploID <- 1:nrow(newhaplo)
newhaplo$Var1 <- as.character(newhaplo$Var1)

#~~~~~~~~
plothaplo <- NULL

for(i in 1:nrow(newhaplo)){
  x <- data.frame(Allele = strsplit(newhaplo$Var1[i], split = "")[[1]],
                  SNP.Name = imputed.snplist$SNP,
                  HaploID = newhaplo$HaploID[i])
  x$Order <- 1:nrow(x)
  plothaplo <- rbind(plothaplo, x)
  rm(x)
}

plothaplo$Allele <- toupper(plothaplo$Allele)

head(plothaplo)
head(imputed.snplist)

testfreq <- imputed.snplist[,1:4]
head(testfreq)
testfreq$refallele <- ifelse(testfreq$Freq1 <= 0.5, testfreq$Al2, testfreq$Al1)

plothaplo$FillColour <- NA
for(i in 1:nrow(plothaplo)){
  plothaplo$FillColour[i] <- ifelse(plothaplo$Allele[i] == testfreq$refallele[which(testfreq$SNP == plothaplo$SNP.Name[i])], "A", "B")
}


importantsnps <- rbind(subset(chr6results, Array == "snp50"),
                       subset(chr6results, SNP.Name == chr6results$SNP.Name[which(chr6results$Pc1df == min(chr6results$Pc1df))]))
importantsnps <- join(importantsnps, unique(select(plothaplo, SNP.Name, Order)))

ggplot() +
  geom_tile(data = plothaplo, aes(Order, HaploID, fill = FillColour), colour = "black") +
  scale_fill_brewer(palette = "Set1") +
  geom_point(data = importantsnps, aes(x = Order, y = 0, colour = Array))



