# Imputation of genotypes using MACH
# Author: Susan Johnston

# Requires plink in PATH, mach1 in the mach directory


setwd("mach")

library(GenABEL)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

#~~ make dat file


system("cmd", input = "plink --file chr6segmentmach --from s15515.1 --to s74824.1 --no-pheno --sheep --recode --out chr6segmentmachsmall")
system("cmd", input = "plink --file chr6HDsegmentmach --from s15515.1 --no-pheno --to oar3_OAR6_117013262 --sheep --recode --out chr6HDsegmentmachsmall")

#~~ Load pedigree

firstrun <- FALSE
if(firstrun == TRUE){
  
  pedigree <- read.table("../data/pedigree_20130920.txt", header = T)
  names(pedigree)[1] <- "V2"
  pedigree$FATHER[which(is.na(pedigree$FATHER))] <- 0
  pedigree$MOTHER[which(is.na(pedigree$MOTHER))] <- 0
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ 1. Create dataset for SNP50 chip    #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~ Create .dat file
  
  ped50.map <- read.table("chr6segmentmachsmall.map", stringsAsFactors = F)
  ped50.dat <- cbind("M", ped50.map$V2)
  write.table(ped50.dat, "test50.dat", row.names = F, col.names = F, quote = F)
  
  #~~ Amend .ped file
  ped50 <- read.table("chr6segmentmachsmall.ped", stringsAsFactors = F)
  head(ped50)
  
  if(all(ped50$V6 == -9)) ped50 <- ped50[,-6]
  
  ped50.test <- join(data.frame(V2 = ped50$V2), pedigree)
  head(ped50.test)
  ped50.test$FATHER[which(!ped50.test$FATHER %in% ped50.test$V2)] <- 0
  ped50.test$MOTHER[which(!ped50.test$MOTHER %in% ped50.test$V2)] <- 0
  
  ped50.test$MOTHER[which(ped50.test$FATHER == 0 & ped50.test$MOTHER !=0)] <- 0
  ped50.test$FATHER[which(ped50.test$MOTHER == 0 & ped50.test$FATHER !=0)] <- 0
  
  
  all(ped50.test$V2 == ped50$V2)
  ped50$V3 <- ped50.test$FATHER
  ped50$V4 <- ped50.test$MOTHER
  ped50$V1 <- 1
  
  write.table(ped50, "test50.ped", row.names = F, col.names = F, quote = F)
  
  system("cmd", input = "mach1 -d test50.dat -p test50.ped --rounds 20 --states 200 --phase --interim 5 --sample 5 --compact --prefix snp50set")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ 2. Create dataset for HD chip       #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  pedHD.map <- read.table("chr6HDsegmentmachsmall.map", stringsAsFactors = F)
  pedHD.dat <- cbind("M", pedHD.map$V2)
  write.table(pedHD.dat, "testHD.dat", row.names = F, col.names = F, quote = F)
  
  #~~ Amend .ped file
  pedHD <- read.table("chr6HDsegmentmachsmall.ped", stringsAsFactors = F)
  pedHD$V1 <- 1
  head(pedHD)
  
  
  if(all(pedHD$V6 == -9)) pedHD <- pedHD[,-6]
  
  write.table(pedHD, "testHD.ped", row.names = F, col.names = F, quote = F)
  
  system("cmd", input = "mach1 -d testHD.dat -p testHD.ped --rounds 20 --states 200 --phase --interim 5 --sample 5 --compact --prefix snpHDset")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ 3. Create imputation set            #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ped50headers <- paste(rep(ped50.map$V2, each = 2),
                        rep(c("a", "b"), times = nrow(ped50.map)),
                        sep = ".")
  
  pedHDheaders <- paste(rep(pedHD.map$V2, each = 2),
                        rep(c("a", "b"), times = nrow(pedHD.map)),
                        sep = ".")
  
  names(ped50) <- c("Family", "ID", "Father", "Mother", "Sex", ped50headers)
  
  names(pedHD) <- c("Family", "ID", "Father", "Mother", "Sex", pedHDheaders)
  
  pedHD <- select(pedHD, -Father, -Mother, -Sex, -Family)
  pedHD <- pedHD[,-which(names(pedHD) %in% ped50headers)]
  
  
  impute.set <- join(ped50, pedHD)
  head(impute.set)
  
  impute.set <- cbind(impute.set[,1:5], impute.set[,pedHDheaders])
  for(i in 6:ncol(impute.set)) impute.set[which(is.na(impute.set[,i])),i] <- 0
  
  write.table(impute.set, "testImpute.ped", row.names = F, col.names = F, quote = F)
  write.table(pedHD.dat, "testImpute.dat", row.names = F, col.names = F, quote = F)
  write.table(pedHD.dat[,2], "testImpute.snps", row.names = F, col.names = F, quote = F)
  
  
  #system("cmd", input = "mach1 -d testImpute.dat -p testImpute.ped -s testImpute.snps -h snpHDset --crossover snpHDset.erate --errormap snpHDset.rec --greedy --rounds 10 --mle --prefix imputed")
  system("cmd", input = "mach1 -d test50.dat -p test50.ped -s testImpute.snps -h snpHDset.haplos --crossover snpHDset.erate --errormap snpHDset.rec --greedy --rounds 10 --mle --mldetails --prefix imputed")
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Read in imputed genotypes     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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








pr.melt <- melt(imputed.pr, id.vars = c("Link"))
head(pr.melt)
names(pr.melt)[2:3] <- c("SNP.Name", "Dosage")

ggplot(pr.melt, aes(Dosage)) + geom_histogram(binwidth = 0.05)



geno.melt <- melt(imputed, id.vars = c("Link"))
head(geno.melt)
names(geno.melt)[2:3] <- c("SNP.Name", "Genotype")

restab <- join(geno.melt, pr.melt)
head(restab)

restab$id <- gsub("1->", "", restab$Link) 

predictors <- read.table("../results/2_predictors_grm.summ.RR.txt", header = T)

#predictors.f <- read.table("../results/2_predictors_grm.summ.RR.females.txt", header = T)
head(predictors)

restab <- join(restab, predictors)

head(restab)
restab$Genotype <- as.character(restab$Genotype)
restab$Genotype[which(restab$Dosage > 0.25 & restab$Dosage < 0.75)] <- "0/0"
restab$Genotype[which(restab$Dosage > 1.25 & restab$Dosage < 1.75)] <- "0/0"

restab2 <- droplevels(filter(restab, !is.na(giv.BLUP)))

table(restab2$Genotype)
unique(restab2[which(restab2$Genotype == "0/0"),"id"])
restab2 <- subset(restab2, !id %in% unique(restab2[which(restab2$Genotype == "0/0"),"id"]))

ggplot(restab2, aes(Dosage)) + geom_histogram(binwidth = 0.05)
ggplot(restab2, aes(Dosage)) + geom_histogram(binwidth = 0.05) + coord_cartesian(ylim = c(0, 100))



modelres <- data.frame(SNP.Name = unique(restab2$SNP.Name), adjR2 = NA, Fstat = NA, df1 = NA, df2 = NA)
modelres$SNP.Name <- as.character(modelres$SNP.Name)

modeleffects <- NULL

for(i in 1:nrow(modelres)){
  fit1 <- lm(Resid.Mean ~ Genotype, data = droplevels(filter(restab2, SNP.Name ==  modelres$SNP.Name[i])))
  modelres$adjR2[i] <- summary(fit1)$adj.r.squared
  modelres[i, c("Fstat", "df1", "df2")] <- summary(fit1)$fstatistic
  modeleffects <- rbind(modeleffects, data.frame(cbind(SNP.Name =  modelres$SNP.Name[i], summary(fit1)$coefficients)))
  rm(fit1)
}
modelres$P <- pf(modelres$Fstat, modelres$df1, modelres$df2, lower.tail=F)

head(imputed.snplist)
names(imputed.snplist)[1] <- "SNP.Name"

modelres <- join(modelres, imputed.snplist)
head(modelres)
modelres$Order <- 1:nrow(modelres)

ggplot(modelres, aes(-log10(P), Rsq)) + geom_point()
ggplot(modelres, aes(Order, -log10(P))) + geom_point()

modeleffects[which(modeleffects$SNP.Name == modelres$SNP.Name[which(modelres$P == min(modelres$P))]),]

fit1 <- lm(Resid.Mean ~ Genotype, data = droplevels(filter(restab2, SNP.Name ==  modelres$SNP.Name[which(modelres$P == min(modelres$P))])))

write.table(modelres, "ImputationAssociationResults.txt", row.names = F, sep = "\t", quote = F)




#~~ Save genos as a PLINK format file

sheepabel <- load.gwaa.data(phe = "../data/1_GenABEL_hornpheno20150126.txt",
                            gen = "../data/1_GenABEL_sheepabelFullQC20150129.gen")

imputed[1:10,1:10]
HDmap <- read.table("chr6HDsegmentmach.map", stringsAsFactors = F)

maptab <- HDmap[which(HDmap$V2 %in% imputed.snplist$SNP.Name),]
all(maptab$V2 == imputed.snplist$SNP.Name)
maptab <- maptab[,-3]


imputed$Link <- gsub("1->", "", imputed$Link)

for(i in 2:ncol(imputed)) imputed[,i] <- gsub("/", " ", imputed[,i], fixed = T)

imputed <- cbind(Family = 1,
                 imputed[,1],
                 MOTHER = 0, FATHER = 0,
                 SEX = 0, 
                 imputed[,2:ncol(imputed)])

temp <- data.frame(id = imputed[,2])
temp <- join(temp, phdata(sheepabel)[,c("id", "sex")])

imputed$SEX <- temp$sex

write.table(imputed, "imputedPLINK.ped", row.names = F, col.names = F, quote = F)
write.table(maptab, "imputedPLINK.map", row.names = F, col.names = F, quote = F)

convert.snp.ped(pedfile = "imputedPLINK.ped", 
                mapfile = "imputedPLINK.map", 
                outfile = "imputedPLINK.genabel", 
                traits = 0, mapHasHeaderLine = F)

phenofile <- imputed[,c(2, 5)]
names(phenofile) <- c("id", "sex")

write.table(phenofile, "imputedPLINK.phenofile", row.names = F, quote = F)


impabel <- load.gwaa.data(genofile = "imputedPLINK.genabel", phenofile = "imputedPLINK.phenofile")

impmerge <- merge.gwaa.data(sheepabel, impabel)
nsnps(impmerge)
nids(impmerge)

#~~ remove bad ids 

impmerge <- impmerge[unique(restab2$id),]

#~~ remove bad loci
locus.accuracy <- read.table("10foldXvalidationAccuracy.Locus.txt", header = T)
badsnps <- as.character(locus.accuracy$SNP.Name[which(locus.accuracy$ProportionMatch < 0.95)])

impmerge <- impmerge[,which(!snp.names(impmerge) %in% badsnps)]


predtab <- read.table("../results/2_predictors_grm.summ.RR.txt", header = T)
predtab.m <- read.table("../results/2_predictors_grm.summ.RR.males.txt", header = T)
predtab.f <- read.table("../results/2_predictors_grm.summ.RR.females.txt", header = T)

source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")
source("C:/Users/Susan Johnston/Desktop/R Functions/GenABELPlotFunctions.R")

testgw <- add.phdata(impmerge, predtab)
gwaa.resid <- qtscore(Resid.Mean ~ sex, data = testgw)

multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]




#~~ Male sheep
testgw <- add.phdata(impmerge, predtab.m)
gwaa.resid <- qtscore(Resid.Mean, data = testgw)

multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]



#~~ Female sheep
testgw <- add.phdata(impmerge, predtab.f)
gwaa.resid <- qtscore(Resid.Mean, data = testgw)

multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, P1df)[1:10,]





genes <- read.table("C:/Users/Susan Johnston/Desktop/Sheep Genome Oar_v3.1/GeneListv1.78.txt", header = T, sep = "\t", stringsAsFactors = F)

geneplot <- subset(genes, Chromosome == "6" & Start > 115700000)
geneplot <- rbind(geneplot, c(6, NA, 116442905, 116430959, NA, NA, NA, "RNF212", NA, NA))
geneplot$Start <- as.numeric(geneplot$Start)
geneplot$Stop  <- as.numeric(geneplot$Stop )
geneplot$MidPos <- apply(geneplot[,c("Start", "Stop")], 1, mean)

geneplot <- arrange(geneplot, Start)
genemelt <- melt(cbind(id = 1:nrow(geneplot), geneplot[,c("Start", "Stop")]), id.vars = "id")



chr6results <- droplevels(subset(gwaa.results, SNP.Name %in% snp.names(impabel)))
chr6results$Array <- "snp50"
chr6results$Array[grep("oar3_", chr6results$SNP.Name)] <- "snpHD"

Pthresh <- 0.05/22273.61

chr6results$Size <- ifelse(chr6results$Array == "snp50", 2, 1)

write.table(chr6results, "ImputationDataForPlotting.txt", row.names = F, sep = "\t", quote = F)

ggplot() + 
  geom_text(data = geneplot, aes(x = MidPos, y = 9, label = Gene.Name), angle = 270, size = 4) +
  geom_line(data = genemelt, aes(value, 8, group = id), size = 6) +
  geom_point(data = chr6results, aes(Position, -log10(Pc1df), colour = Array, shape = Array), size = 4, alpha = 0.7) +
  scale_colour_manual(values = c("black", "red")) +
  geom_point(data = subset(chr6results, Array == "snp50"), aes(Position, -log10(Pc1df)), size = 5, shape=21, fill = "black") +
  geom_hline(yintercept = -log10(Pthresh),linetype = 2, alpha = 0.6, size = 1) +
  scale_x_continuous(breaks = seq(0, 118000000, 200000), labels = seq(0, 118, 0.2)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x = "Position (MB)", y = "-log10(P)") +
  coord_cartesian(ylim = c(-0.2, 10))


#~~ Make file for haploview
write.table(snpnames(impabel), "snpnamespostQC.txt", row.names = F, col.names = F, quote = F)
system("plink --file chr6HDsegmentmach --no-pheno --extract snpnamespostQC.txt --sheep --recode --out chr6HDsegPostQC")
infofile <- read.table("chr6HDsegPostQC.map")
head(infofile)
write.table(infofile[,c(2, 4)], "chr6HDsegPostQC.info", row.names = F, quote = F, col.names = F)
pedfile <- read.table("chr6HDsegPostQC.ped")
head(pedfile)[,1:10]
pedfile[,6] <- 0
write.table(pedfile, "chr6HDsegPostQC.ped", row.names = F, quote = F, col.names = F)


#~~ Plot without gene info

ggplot() + 
  geom_point(data = chr6results, aes(Position, -log10(Pc1df), colour = Array, shape = Array), size = 2) +
  scale_colour_manual(values = c("black", "red")) +
  geom_point(data = subset(chr6results, Array == "snp50"), aes(Position, -log10(Pc1df)), size = 2, col = "black") +
  geom_hline(yintercept = -log10(Pthresh),linetype = 2, alpha = 0.6, size = 0.5) +
  scale_x_continuous(breaks = seq(0, 118000000, 200000), labels = seq(0, 118, 0.2)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, vjust = 0),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12, vjust = 0.7),
        axis.title.y = element_text (size = 12, angle = 90, vjust = 0.7),
        axis.title.x = element_text (size = 12, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x = "Position (MB)", y = "-log10(P)") +
  coord_cartesian(xlim = c(115700000, 117031472))


#~~ LD plot

sheepHD <- load.gwaa.data(phe = "../../Soay Sheep HD SNP Chip/data_from_WTCRF/E10905_SheepHD_Plates1-2_310114 QC1/20140214_SheepHD_QC1_GenABELpheno.txt",
                          gen = "../../Soay Sheep HD SNP Chip/data_from_WTCRF/E10905_SheepHD_Plates1-2_310114 QC1/20140214_SheepHD_QC1_GenABEL.txt")

subHD <- sheepHD[,which(chromosome(sheepHD) == 6 & map(sheepHD) > 115716891)]
nsnps(subHD)
nsnps(impabel)
subHD <- subHD[,snpnames(impabel)]


genoHD    <- data.frame(as.character.gwaa.data(subHD))
genoHDimp <- data.frame(as.character.gwaa.data(impabel[which(idnames(impabel) %in% unique(restab2$id)),]))



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


print(p, vp=viewport(angle=-90))
ggplot(ldtab, aes(Allele1, Allele2, fill = HDcor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

print(p, vp=viewport(angle=-45))


library(LDheatmap)

colfunc <- colorRampPalette(c("red", "white"))
colfunc(20)

plot(rep(1,20),col=colfunc(20),pch=19,cex=3)

str(p)

impabel.genotab <- as.genotype.gwaa.data(impabel[,which(snp.names(impabel) %in% chr6results$SNP.Name)])
p <- LDheatmap::LDheatmap(impabel.genotab, 
                     genetic.distances = map(impabel[,which(snp.names(impabel) %in% chr6results$SNP.Name)]),
                     flip = T,
                     color = colfunc(20))

save(p, file = "ldplot.Rdata")
beepr::beep()

load("ldplot.Rdata")
plot(p)
str(p)

test <- data.frame(p$LDmatrix)
head(test)
test[1:5,1:5]
test$SNP.Name <- row.names(test)

test <- melt(test, id.vars = "SNP.Name")
test <- subset(test, !is.na(value))

test$variable <- as.character(test$variable)

test <- test[unique(sort(c(which(test$SNP.Name == "s74824.1"), which(test$variable == "s74824.1")))),]
test$SNP.Name[which(test$SNP.Name == "s74824.1")] <- test$variable[which(test$SNP.Name == "s74824.1")]
test <- select(test, -variable)
chr6results <- join(chr6results, test)



ggplot() + 
  geom_point(data = chr6results, aes(Position, -log10(Pc1df), colour = Array, shape = Array), size = 2) +
  scale_colour_manual(values = c("black", "red")) +
  geom_point(data = subset(chr6results, Array == "snp50"), aes(Position, -log10(Pc1df)), size = 2, col = "black") +
  geom_point(data = chr6results, aes(x = Position, y = 8, col = value), size = 2) +
  geom_hline(yintercept = -log10(Pthresh),linetype = 2, alpha = 0.6, size = 0.5) +
  scale_x_continuous(breaks = seq(0, 118000000, 200000), labels = seq(0, 118, 0.2)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, vjust = 0),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12, vjust = 0.7),
        axis.title.y = element_text (size = 12, angle = 90, vjust = 0.7),
        axis.title.x = element_text (size = 12, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x = "Position (MB)", y = "-log10(P)") +
  coord_cartesian(xlim = c(115700000, 117031472))


ggplot(subset(chr6results, !is.na(value)), aes(x = Position, y = 8.2)) +
  geom_point(aes(col = value), shape = 15) +
  scale_colour_gradient(low = "white", high = "red") +
  theme_bw()

write.table(chr6results, "test.txt", row.names = F, sep = "\t", quote = F)

head(chr6results)
chr6results <- arrange(chr6results, Position)
chr6results$Sex <- "Female"
chr6results <- select(chr6results, Sex, SNP.Name, Chromosome, Position, A1, A2, N, effAB, effBB, Pc1df, Array, value)
                      
HDstats <- summary.snp.data(gtdata(sheepHD[,which(chromosome(sheepHD) == 6)]))
head(HDstats)
HDstats$SNP.Name <- row.names(HDstats)
HDstats <- select(HDstats, SNP.Name, CallRate, Q.2)

chr6results <- join(chr6results, HDstats)
head(chr6results)
dput(names(chr6results))

names(chr6results) <- c("Sex", "SNP.Name", "Chromosome", "Position.bp", "Allele1", "Allele2", "N", 
                        "EffectAB", "EffectBB", "P.Corrected", "Array", "r2.with.s74824.1", "Call.Rate", "Minor.Allele.Freq")

write.table(chr6results, "C:/Users/Susan Johnston/Desktop/Recombination Rate Manuscript/SupplementaryInformation/ImputationSignificance.csv",
            row.names = F, quote = F, sep = ", ")                      
                      
# test <- dprfast(impabel)
# row.names(test)
# test2 <- test[lower.tri(test)]
# test2 <- melt(test)
# 
# source("C:/Users/Susan Johnston/Desktop/R Functions/recoderFunc.R")
# test2$Allele1 <- as.numeric(recoderFunc(test2$X1, snp.names(impabel), 1:length(snp.names(impabel))))
# test2$Allele2 <- as.numeric(recoderFunc(test2$X2, snp.names(impabel), 1:length(snp.names(impabel))))
# 
# test2 <- (test2[which(test2$Allele2 > test2$Allele1),])
# 
# 
# 
# ggplot(test2, aes(Allele2, Allele1, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red")


#~~ Redo for the different genotypes

genoHD.AA <- subset(data.frame(as.character.gwaa.data(subHD)), s74824.1 == "A/A")
genoHD.GA <- subset(data.frame(as.character.gwaa.data(subHD)), s74824.1 == "G/A")
genoHD.GG <- subset(data.frame(as.character.gwaa.data(subHD)), s74824.1 == "G/G")


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

i = 80









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






head(imputed.haplos)

table(imputed.haplos$V3)
imputed.haplos$V3 <- toupper(imputed.haplos$V3)

imputed.haplos
imputed.haplos$V4 <- unlist(sapply(imputed.haplos$V3, function(x) substr(as.character(x), 45, 126)[[1]]))

V4tab <- data.frame(table(imputed.haplos$V4))
head(V4tab)
V4tab$HaploGroup <- letters[1:nrow(V4tab)]
head(imputed.haplos)
names(V4tab)[1] <- "V4"

imputed.haplos <- join(imputed.haplos, V4tab)

#~~~~~~~~
plothaplo <- NULL

for(i in 1:nrow(imputed.haplos)){
  x <- data.frame(Allele = strsplit(as.character(imputed.haplos$V4[i]), split = "")[[1]],
                  SNP.Name = imputed.snplist$SNP[45:126],
                  HaploID = i,
                  HaploGroup = imputed.haplos$HaploGroup[i])
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

#~~ most common haplotype

testfreq$CommonAllele <- NA
testfreq$CommonAllele[45:126] <-   strsplit(as.character(V4tab$V4[which(V4tab$Freq == max(V4tab$Freq))]), split = "")[[1]]


# plothaplo$FillColour <- NA
# for(i in 1:nrow(plothaplo)){
#   plothaplo$FillColour[i] <- ifelse(plothaplo$Allele[i] == testfreq$refallele[which(testfreq$SNP == plothaplo$SNP.Name[i])], "A", "B")
# }

plothaplo$FillColour <- NA
for(i in 1:nrow(plothaplo)){
  plothaplo$FillColour[i] <- ifelse(plothaplo$Allele[i] == testfreq$CommonAllele[which(testfreq$SNP == plothaplo$SNP.Name[i])], "A", "B")
}


importantsnps <- rbind(subset(chr6results, Array == "snp50"),
                       subset(chr6results, SNP.Name == chr6results$SNP.Name[which(chr6results$Pc1df == min(chr6results$Pc1df))]))
importantsnps <- join(importantsnps, unique(select(plothaplo, SNP.Name, Order)))


plothaplo$HaploGroup <- factor(plothaplo$HaploGroup, levels = rev(c("l", "n", "c", "j", "g", "f", "b", "m", "h", "i", "e", "k", "d", "a")))
plothaplo <- arrange(plothaplo, HaploGroup)
plothaplo$NewOrder <- rep(1:376, each =  82)


ggplot() +
  geom_tile(data = plothaplo, aes(Order, NewOrder, fill = FillColour), colour = "black") +
  scale_fill_brewer(palette = "Set1") +
  geom_point(data = importantsnps, aes(x = Order, y = 0, colour = Array))

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
  



