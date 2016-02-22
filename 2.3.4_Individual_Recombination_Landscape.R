#
#
# Characterisation of the recombination landscape
# Individual measures
#
# Susan Johnston
#
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Binsize we will work with

binsize <- 1000000

#~~ load functions and libraries

library(GenABEL)
library(beepr)
library(ggplot2)
library(reshape)
library(plyr)
library(dplyr)
library(beepr)
library(asreml)

options(stringsAsFactors = F)

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in data

# rectab

rectab <- read.table(paste0("results/2_IndivRR_FullCleanPostSim_", AnalysisSuffix, ".txt"), header = T)
recsumm <- read.table(paste0("results/2_TotalSummIndivRR_FullCleanPostSim_", AnalysisSuffix, ".txt"), header = T)


# Map data 

maptab <- read.table(paste0("results/2_Merged_map_", AnalysisSuffix, ".txt"), header = T)
tail(maptab)
maptab$Cumu <- cumsum(c(maptab$GenomicPosition[1], as.numeric(maptab$GenomicDiff)))[1:(nrow(maptab))]
maxvals <- read.table(paste0("results/2_MaxVals_", AnalysisSuffix, ".txt"), header = T)

# defopts

defopts <- theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 strip.background = element_blank())


#~~ Genomic data

sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")

imputed.genos <- load.gwaa.data(phe = "mach/imputedPLINK.phenofile",
                                gen = "mach/imputedPLINK.genabel")


#~~ Crossover prs population wide

probtab <- read.table(paste0("results/2_PopWide_CrossoverPr_bin", binsize/1000, "K_Full_", AnalysisSuffix, ".txt"), header = T)
probtab$BinID <- paste(probtab$Chr, probtab$BinID, sep = "_")


#~~ asreml stuff

pedigree <- read.table("data/pedigree_20130920.txt", header = T)

names(pedigree)[1] <- "RRID"
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])
dim(pedigree)

ainv <- asreml.Ainverse(pedigree)$ginv


load("gcta/Autosomes_GRM.Rdata")

source("r/makeGRM.R")
source("r/ASReml.EstEffects.R")
source("r/ASReml.ExtractPredictors.R")
source("r/pin.R")
source("r/multiplot.R")
source("r/GenABELPlotFunctions.R")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Create function to determine recomb probabilities #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

extractCrossBin <- function(string, map.vec, binsize, chr, identifier = NULL){
  
  df <- data.frame(Phase = strsplit(string, split = "")[[1]],
                   Position = map.vec)
  df$Phase[df$Phase == "i"] <- 1
  df$Phase[df$Phase == "o"] <- 0
  df$Phase[df$Phase == ":"] <- "-"
  df$Phase[df$Phase == "c"] <- "-"
  df <- filter(df, Phase != "-")
  
  
  df$Temp <- c(-9, df$Phase[-nrow(df)])
  xovertab <- data.frame(Start = df$Position[which(df$Phase != df$Temp)[-1] - 1],
                         Stop = df$Position[which(df$Phase != df$Temp)[-1]])
  
  
  # Make frame for the bins which are spanned.
  
  if(nrow(xovertab) > 0){
    
    bintab <- NULL
    
    for(i in 1:nrow(xovertab)){
      
      tempspan <- as.numeric(as.character(cut(unlist(xovertab[i,]),  seq(1, max(xovertab[i,]) + binsize, binsize), labels = seq(1, max(xovertab[i,]), binsize))))
      tempbintab <- data.frame(Bin = seq(tempspan[1], tempspan[2], binsize))
      tempbintab$BinInitProb <- NA
      
      
      if(nrow(tempbintab) == 1) tempbintab$BinInitProb <- 1
      
      
      if(nrow(tempbintab) > 1){
        tempbintab$BinInitProb[1] <- ((tempbintab$Bin[1] + binsize - 1) - xovertab[i,1])/binsize
        tempbintab$BinInitProb[nrow(tempbintab)] <- (xovertab[i,2] - tempbintab$Bin[nrow(tempbintab)])/binsize
        tempbintab$BinInitProb[which(is.na(tempbintab$BinInitProb))] <- 1
        
      }
      
      totalxoverdist <- sum(tempbintab$BinInitProb)
      tempbintab$BinProb <- tempbintab$BinInitProb/sum(tempbintab$BinInitProb)
      if(!is.null(identifier)) tempbintab$identifier <- identifier
      tempbintab$Chr <- chr
      
      bintab <- rbind(bintab, tempbintab)
      
      rm(tempspan, tempbintab, totalxoverdist)
    }
    
    bintab
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Calculate individual recombination probabilities #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

RunMasterbin == FALSE
if(RunMasterbin == TRUE){

masterbin <- list()

system.time({
  for(i in 1:nrow(rectab)){
    
    if(i %in% seq(1, nrow(rectab), 1000)) print(paste("Analysing row", i, "of", nrow(rectab)))
    
    masterbin[[i]] <- extractCrossBin(rectab$data.v2[i],
                                      map.vec = maptab$GenomicPosition[which(maptab$Chr == rectab$Chr[i])],
                                      binsize = binsize,
                                      chr = rectab$Chr[i],
                                      identifier = i)
    
  }
})

system.time(masterbin2 <- data.table::rbindlist(masterbin))

head(masterbin2)

#~~ Merge with IDs

masterbin2 <- join(masterbin2,
                   data.frame(UniqueID = rectab$UniqueID,
                              UniqueID2 = rectab$UniqueID2,
                              RRID = rectab$RRID,
                              RRID.SEX = rectab$RRID.SEX,
                              identifier = 1:nrow(rectab)))

#~~ Give Bins more information

masterbin2 <- data.frame(masterbin2)
masterbin2 <- subset(masterbin2,  select = c(Bin, BinProb, Chr, UniqueID, UniqueID2, RRID, RRID.SEX))



write.table(masterbin2,
            paste0("results/2_Indiv_CrossoverPr_bin", binsize/1000, "K_", AnalysisSuffix, ".txt"),
            row.names = F, sep = "\t", quote = F)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Calculate telomeric crossover count              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


masterbin2 <- read.table(paste0("results/2_Indiv_CrossoverPr_bin", binsize/1000, "K_", AnalysisSuffix, ".txt"), header = T)

idtemplate <- unique(dplyr::select(rectab, UniqueID2, RRID, RRID.SEX))

tempbin <- join(masterbin2, idtemplate)

tempbin <- join(tempbin, probtab[,c("Bin", "Chr", "BinID", "DistanceToTelo", "PC.GC", "SNP.Count", "TelomereClass")])
head(tempbin)

rnfgeno <- data.frame(rnf212 = as.character.gwaa.data(imputed.genos[,"oar3_OAR6_116402578"]))
rnfgeno$RRID <- row.names(rnfgeno)
names(rnfgeno)[1] <- "rnf212"

tempbin <- join(tempbin, rnfgeno)
idtemplate <- join(idtemplate, rnfgeno)

table(idtemplate$rnf212)
sex.geno.count <- data.frame(table(idtemplate$rnf212, idtemplate$RRID.SEX))
sex.geno.count$Name <- paste(sex.geno.count$Var1, sex.geno.count$Var2, sep = ".")
sex.geno.count$Name <- gsub("/", ".", sex.geno.count$Name)


genobin <- data.frame(tapply(tempbin$BinProb, list(tempbin$BinID, tempbin$rnf212, tempbin$RRID.SEX), sum, na.rm = T))
head(genobin)

for(i in 1:ncol(genobin)){
  genobin[,i] <- genobin[,i]/sex.geno.count$Freq[which(sex.geno.count$Name == names(genobin[i]))]
}

genobin$BinID <- row.names(genobin)

genobin <- join(genobin, probtab)
head(genobin)

genobin <- select(genobin, -MaxBin, -DistanceToEnd, -Count.A, -Count.C, -Count.G, -Count.T, -Count.ACGT)

genobin.sexgeno <- melt(genobin, id.vars = c("BinID", "Bin", "Chr", "Sex.averaged.Pr", "Male.Pr", 
                                             "Female.Pr", "DistanceToTelo", "MaleToFemale", "SNP.Count", "Mean.MAF", 
                                             "Mean.Inf.Mei", "PC.GC", "AdjustedDistance", "TelomereClass"))

head(genobin.sexgeno)

genobin$AA.Ratio <- genobin$A.A.Male/genobin$A.A.Female
genobin$GA.Ratio <- genobin$G.A.Male/genobin$G.A.Female
genobin$GG.Ratio <- genobin$G.G.Male/genobin$G.G.Female

genobin.ratiomelt <- genobin.sexgeno <- melt(genobin, id.vars = c("BinID", "Bin", "Chr", "Sex.averaged.Pr", "Male.Pr", 
                                                                  "Female.Pr", "DistanceToTelo", "MaleToFemale", "SNP.Count", "Mean.MAF", 
                                                                  "Mean.Inf.Mei", "PC.GC", "AdjustedDistance", "TelomereClass"))

head(genobin.ratiomelt)
table(genobin.ratiomelt$variable)
genobin.ratiomelt <- droplevels(genobin.ratiomelt[grep("Ratio", genobin.ratiomelt$variable),])


ggplot(subset(genobin.sexgeno, DistanceToTelo < 61), aes(DistanceToTelo,  value, col = variable)) +
  geom_point(alpha = 0.2) +
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
        legend.background = element_blank())

genobin.ratiomelt[which(genobin.ratiomelt$value > 50),]


save(genobin.sexgeno, genobin.ratiomelt, file = "results/2_DataForGenoSpecificRecombLandscapev2.Rdata")

multiplot(
ggplot(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Female", "G.A.Female", "G.G.Female")),
       aes(DistanceToTelo,  value, col = variable)) +
  geom_point(alpha = 0.2) +
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
        legend.background = element_blank())
,
ggplot(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Male", "G.A.Male", "G.G.Male")),
       aes(DistanceToTelo,  value, col = variable)) +
  geom_point(alpha = 0.2) +
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
        legend.background = element_blank())
,
ggplot(filter(genobin.ratiomelt, DistanceToTelo < 61, value < 50), aes(DistanceToTelo,  value, col = variable)) +
  geom_point(alpha = 0.2) +
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
        legend.background = element_blank())

, cols = 3)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Calculate telomeric crossover count              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


masterbin2 <- read.table(paste0("results/2_Indiv_CrossoverPr_bin", binsize/1000, "K_", AnalysisSuffix, ".txt"), header = T)
idtemplate <- unique(dplyr::select(rectab, UniqueID2, RRID, RRID.SEX))


tempbin <- join(masterbin2, probtab[,c("Bin", "Chr", "DistanceToTelo")])
head(tempbin)


telobin    <- filter(tempbin, DistanceToTelo <= 10)
nontelobin <- filter(tempbin, DistanceToTelo >  11)

telocount <- data.frame(NonTelXover = tapply(telobin$BinProb, telobin$UniqueID2, sum))
telocount$UniqueID2 <- row.names(telocount)

nontelocount <- data.frame(NonNonTelXover = tapply(nontelobin$BinProb, nontelobin$UniqueID2, sum))
nontelocount$UniqueID2 <- row.names(nontelocount)

telorecsumm <- join(idtemplate, telocount)
telorecsumm <- join(telorecsumm, nontelocount)


telorecsumm$NonTelXover[which(is.na(telorecsumm$NonTelXover))] <- 0
telorecsumm$NonNonTelXover [which(is.na(telorecsumm$NonNonTelXover))] <- 0

telorecsumm <- join(telorecsumm, recsumm[,c("UniqueID2", "TotalRecombCount", "RRID.Fhat3")])


ggplot(telorecsumm, aes(NonTelXover))    + geom_histogram(binwidth = 1) + facet_wrap(~RRID.SEX, ncol = 1)
ggplot(telorecsumm, aes(NonNonTelXover)) + geom_histogram(binwidth = 1) + facet_wrap(~RRID.SEX, ncol = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Is there a heritable basis to this?              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


telorecsumm$RRID <- as.factor(telorecsumm$RRID)
telorecsumm$RRID.SEX <- as.factor(telorecsumm$RRID.SEX)

attributes(telorecsumm$NonTelXover) <- NULL
attributes(telorecsumm$NonNonTelXover) <- NULL

telorecsumm.m <- droplevels(filter(telorecsumm, RRID.SEX == "Male"))
telorecsumm.f <- droplevels(filter(telorecsumm, RRID.SEX == "Female"))

# ALL IDS

grminv <- makeGRM(grm.auto, ids.auto, telorecsumm$RRID)

telo.summ.RR <- asreml(fixed = NonTelXover ~ RRID.SEX + RRID.Fhat3,
                      random = ~ giv(RRID) + ide(RRID),
                      data = telorecsumm,
                      ginverse =  list(RRID = grminv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

summary(telo.summ.RR, all = T)$coef.fixed
ASReml.EstEffects(telo.summ.RR)
wald.asreml(telo.summ.RR)



# MALES ONLY
grminv <- makeGRM(grm.auto, ids.auto, telorecsumm.m$RRID)

telo.summ.RR.males <- asreml(fixed = TelXover ~ RRID.Fhat3,
                            random = ~ giv(RRID) + ide(RRID),
                            data = telorecsumm.m,
                            ginverse =  list(RRID = grminv),
                            na.method.X = "omit", na.omit.Y = "na.omit",
                            workspace = 500e+6, pworkspace = 500e+6)


summary(telo.summ.RR.males, all = T)$coef.fixed
wald.asreml(telo.summ.RR.males)
ASReml.EstEffects(telo.summ.RR.males)


# FEMALES ONLY
grminv <- makeGRM(grm.auto, ids.auto, telorecsumm.f$RRID)

telo.summ.RR.females <- asreml(fixed = TelXover ~ RRID.Fhat3,
                              random = ~ giv(RRID) + ide(RRID),
                              data = telorecsumm.f,
                              ginverse =  list(RRID = grminv),
                              na.method.X = "omit", na.omit.Y = "na.omit",
                              workspace = 500e+6, pworkspace = 500e+6)

summary(telo.summ.RR.females, all = T)$coef.fixed
wald.asreml(telo.summ.RR.females)
ASReml.EstEffects(telo.summ.RR.females)


# ALL IDS

grminv <- makeGRM(grm.auto, ids.auto, telorecsumm$RRID)

nontelo.summ.RR <- asreml(fixed = NonTelXover ~ RRID.SEX + RRID.Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = telorecsumm,
                       ginverse =  list(RRID = grminv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(nontelo.summ.RR, all = T)$coef.fixed
ASReml.EstEffects(nontelo.summ.RR)
wald.asreml(nontelo.summ.RR)


# MALES ONLY
grminv <- makeGRM(grm.auto, ids.auto, telorecsumm.m$RRID)

nontelo.summ.RR.males <- asreml(fixed = NonTelXover ~ RRID.Fhat3,
                             random = ~ giv(RRID) + ide(RRID),
                             data = telorecsumm.m,
                             ginverse =  list(RRID = grminv),
                             na.method.X = "omit", na.omit.Y = "na.omit",
                             workspace = 500e+6, pworkspace = 500e+6)



summary(nontelo.summ.RR.males, all = T)$coef.fixed
wald.asreml(nontelo.summ.RR.males)
ASReml.EstEffects(nontelo.summ.RR.males)


# FEMALES ONLY
grminv <- makeGRM(grm.auto, ids.auto, telorecsumm.f$RRID)

nontelo.summ.RR.females <- asreml(fixed = NonTelXover ~ RRID.Fhat3,
                               random = ~ giv(RRID) + ide(RRID),
                               data = telorecsumm.f,
                               ginverse =  list(RRID = grminv),
                               na.method.X = "omit", na.omit.Y = "na.omit",
                               workspace = 500e+6, pworkspace = 500e+6)


summary(nontelo.summ.RR.females, all = T)$coef.fixed
wald.asreml(nontelo.summ.RR.females)
ASReml.EstEffects(nontelo.summ.RR.females)


rm(list = grep("\\.x", ls(), value = T))
rm(list = grep("\\.wogiv", ls(), value = T))
gc()

ls()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. GWAS                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Pthresh <- 0.05/22273.61

teloRRvals <- extractPredictors(telo.summ.RR, telorecsumm$RRID, idvar = "RRID")
names(teloRRvals)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, teloRRvals)
gwaa.resid <- qtscore(Resid.Mean ~ sex, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]
FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh)


# Males Only:

teloRRvals.male <- extractPredictors(telo.summ.RR.males, telorecsumm.m$RRID, idvar = "RRID")
names(teloRRvals.male)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, teloRRvals.male)
gwaa.resid <- qtscore(Resid.Mean, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))
gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]


# Females Only

teloRRvals.female <- extractPredictors(telo.summ.RR.females, telorecsumm.f$RRID, idvar = "RRID")
names(teloRRvals.female)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, teloRRvals.female)
gwaa.resid <- qtscore(Resid.Mean, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]

# all, nontelo

nonteloRRvals <- extractPredictors(nontelo.summ.RR, telorecsumm$RRID, idvar = "RRID")
names(nonteloRRvals)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, nonteloRRvals)
gwaa.resid <- qtscore(Resid.Mean ~ sex, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]
FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh)


# Males Only:

nonteloRRvals.male <- extractPredictors(nontelo.summ.RR.males, telorecsumm.m$RRID, idvar = "RRID")
names(nonteloRRvals.male)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, nonteloRRvals.male)
gwaa.resid <- qtscore(Resid.Mean, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))
gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]


# Females Only

nonteloRRvals.female <- extractPredictors(nontelo.summ.RR.females, telorecsumm.f$RRID, idvar = "RRID")
names(nonteloRRvals.female)[1] <- "id"

sheepabel2 <- add.phdata(sheepabel, nonteloRRvals.female)
gwaa.resid <- qtscore(Resid.Mean, data = sheepabel2)
multiplot(FullGwasPlot(CumuPos(gwaa.resid), corrected = T,bonf = Pthresh),
          FullPpPlot(gwaa.resid, CumuPos(gwaa.resid), corrected = T), cols = 2, layout = matrix(c(1, 1, 2), nrow=1, byrow=T))

gwaa.results <- results(gwaa.resid)
gwaa.results$SNP.Name <- row.names(gwaa.results)
arrange(gwaa.results, Pc1df)[1:10,]



# 
# rnfgeno <- data.frame(rnf212 = as.character.gwaa.data(sheepabel[,"s74824.1"]))
# rnfgeno$RRID <- row.names(rnfgeno)
# names(rnfgeno)[1] <- "rnf212"
# 
# telorecsumm.f <- join(telorecsumm.f, rnfgeno)
# telorecsumm.f$rnf212 <- as.factor(telorecsumm.f$rnf212)
# 
# telo.summ.RR.females.rnf <- asreml(fixed = TelXover ~ RRID.Fhat3 + rnf212,
#                                random = ~ giv(RRID) + ide(RRID),
#                                data = telorecsumm.f,
#                                ginverse =  list(RRID = grminv),
#                                na.method.X = "omit", na.omit.Y = "na.omit",
#                                workspace = 500e+6, pworkspace = 500e+6)
# 
# 
# summary(telo.summ.RR.females.rnf, all = T)$coef.fixed
# wald.asreml(telo.summ.RR.females.rnf)
# ASReml.EstEffects(telo.summ.RR.females.rnf)
# 
# nontelo.summ.RR.females.rnf <- asreml(fixed = NonTelXover ~ RRID.Fhat3 + rnf212,
#                                    random = ~ giv(RRID) + ide(RRID),
#                                    data = telorecsumm.f,
#                                    ginverse =  list(RRID = grminv),
#                                    na.method.X = "omit", na.omit.Y = "na.omit",
#                                    workspace = 500e+6, pworkspace = 500e+6)
# 
# 
# summary(nontelo.summ.RR.females.rnf, all = T)$coef.fixed
# wald.asreml(nontelo.summ.RR.females.rnf)
# ASReml.EstEffects(nontelo.summ.RR.females.rnf)
