#
# Results and Figures for Manuscript
#


#+ echo = F, warning = F, message = F
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load data and libraries        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(GenABEL)
library(reshape)
library(plyr)
library(xtable)
library(dplyr)

source("r/multiplot.R")
source("r/recoderFunc.R")


sheepabel <- load.gwaa.data(phe = "data/1_GenABEL_hornpheno20150126.txt",
                            gen = "data/1_GenABEL_sheepabelFullQC20150129.gen")

# recsumm - recombination counts for individuals

recsumm <- read.table("results/2_TotalSummIndivRR_FullCleanPostSim_g.txt", header = T)

# maptab - genomic positions and linkage maps

maptab <- read.table("results/2_Merged_map_g.txt", header = T)
maptab$Cumu <- cumsum(c(maptab$GenomicPosition[1], as.numeric(maptab$GenomicDiff)))[1:(nrow(maptab))]

maptab.sex <- melt(maptab[,c("Chr", "GenomicPosition", "cMPosition.Female", "cMPosition.Male")], id.vars = c("Chr", "GenomicPosition"))
maptab.sex$variable <- gsub("cMPosition.", "", maptab.sex$variable)
maptab.sex$Chr2 <- maptab.sex$Chr
maptab.sex$Chr2 <- gsub("27", "X", maptab.sex$Chr2)

maptab$Cumu2 <- maptab$Cumu + (25000000 * maptab$Chr)


# maxvals - chromosome-wide map lengths

maxvals <- read.table("results/2_MaxVals_g.txt", header = T)
maxvals$Chr2 <- maxvals$Chr
maxvals$Chr2 <- gsub("27", "X", maxvals$Chr2)


sum(round(maxvals$maxGenome[-27]/48e6))

# save(rnf.cistrans.fixef, rnf.cistrans.ranef, file = "results/2_rnf_cistransmodels.Rdata")
# write.table(all.ranef, "results/2_GRMvPed_ModelTableForLaTeX.txt", row.names = F, quote = F, sep = " & ")
# save(rnf.cistrans.fixef, rnf.cistrans.ranef, file = "results/2_rnf_cistransmodels.Rdata")
# fixef.rnf2$Genotype <- rep(c("AA", "AG", "GG"), 3)

#-

head(maptab)
temp <- data.frame(min = tapply(maptab$GenomicPosition, maptab$Chr, min),
                   max = tapply(maptab$GenomicPosition, maptab$Chr, max))
temp <- temp[1:26,]
temp$Diff <- temp$max - temp$min
colSums(temp)



ggplot(recsumm, aes(RRID.SEX, TotalRecombCount, fill = RRID.SEX)) +
  geom_boxplot(col = "grey10", binwidth = 2, notch = T, width = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        plot.title = element_text(size = 18, hjust = 0, vjust = -4),
        #strip.background = element_blank(),
        legend.position = "none")+
  labs(x = "SEX", y = "Total Crossover Count")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 1. Genetic dataset and recombination landscape       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Number of SNPs passing quality control
nsnps(sheepabel)

#' Number of SNPs with known position
nsnps(sheepabel[,chromosome(sheepabel) != 0])

#' Number of IDs and unique measurements
nrow(recsumm); length(unique(recsumm$RRID))
nrow(filter(recsumm, RRID.SEX == "Male")); length(unique(filter(recsumm, RRID.SEX == "Male")$RRID))
nrow(filter(recsumm, RRID.SEX == "Female")); length(unique(filter(recsumm, RRID.SEX == "Female")$RRID))

#' Number of crossovers

sum(recsumm$TotalRecombCount)
sum(filter(recsumm, RRID.SEX == "Male")$TotalRecombCount)
sum(filter(recsumm, RRID.SEX == "Female")$TotalRecombCount)

#' Sex-averaged linkage map lengths

sum(maxvals[which(maxvals$Chr != 27),]$maxcM)
sum(maxvals[which(maxvals$Chr != 27),]$maxcM.male)
sum(maxvals[which(maxvals$Chr != 27),]$maxcM.female)

#' Male to female ratio

sum(maxvals[which(maxvals$Chr != 27),]$maxcM.male)/sum(maxvals[which(maxvals$Chr != 27),]$maxcM.female)

#' Broad scale recombination variation

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

#' Broad scale recombination variation by sex

maptab.sex$Chr2 <- maptab.sex$Chr
maptab.sex$Chr2[which(maptab.sex$Chr2 == 27)] <- "X"
maptab.sex$Chr2 <- factor(maptab.sex$Chr2, levels = c(1:26, "X")) 


plot1 <- ggplot(maptab.sex, aes(GenomicPosition/1e6, value, col = variable)) + 
  geom_line() +
  scale_colour_manual(values = c("black", "darkgrey")) +
  facet_wrap(~Chr2, scale = "free") +
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

pdf("PerChrLinkageMaps.pdf", width = 12, height = 10, useDingbats = F) 
plot1
dev.off()




#' Relationships between chromosome lengths etc.


fit.a <- summary(lm(maxcM ~ maxGenome, data = subset(maxvals, Chr != 27)))
fit.b <- summary(lm(I(maxGenome/maxcM) ~ I(1/maxcM), data = subset(maxvals, Chr != 27)))
fit.c <- summary(lm(maxcM.female ~ maxcM.male, data = subset(maxvals, Chr != 27)))

fit.a$adj.r.squared
fit.a$coefficients[2,4]

fit.b$adj.r.squared
fit.b$coefficients[2,4]

fit.c$adj.r.squared
fit.c$coefficients[2,4]

newpath <- ""

tiff(paste0(newpath, "Fig1.tiff"), width = 7.4, height = 2.5, units = "in", res = 300) 

multiplot(
  ggplot(subset(maxvals, Chr != 27), aes(maxGenome/1e6, maxcM)) + 
    stat_smooth(method = "lm") +
    geom_text(aes(label = Chr2), size = 3) +
    annotate("text", x = maxvals$maxGenome[27]/1e6, y = maxvals$maxcM[27], label = "X", size = 3) +
    labs(x = "Chromosome Length (Mb)", y = "Linkage Map Length (cM)") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 8, vjust = 0),
          axis.text.y  = element_text (size = 8, hjust = 1.3),
          strip.text.x = element_text (size = 8, vjust = 0.7),
          axis.title.y = element_text (size = 8, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 8, vjust = 0.2),
          plot.title = element_text(size = 10, hjust = 0),
          strip.background = element_blank()) +
    ggtitle("A") + 
    scale_x_continuous(breaks = seq(0, 300, 50))
  ,
  
  ggplot(subset(maxvals, Chr != 27), aes(maxGenome/1e6, maxcM/maxGenome*1e6)) + 
    stat_smooth(method = "lm", formula = y ~ I(1/x)) +
    geom_text(aes(label = Chr2), size = 3) +
    annotate("text", x = maxvals$maxGenome[27]/1e6, y = maxvals$maxcM[27]/maxvals$maxGenome[27]*1e6, label = "X", size = 3) +
    coord_cartesian(ylim = c(0.7, 1.9)) +
    labs(x = "Chromosome Length (Mb)", y = "Recombination Rate (cM/Mb)") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 8, vjust = 0),
          axis.text.y  = element_text (size = 8, hjust = 1.3),
          strip.text.x = element_text (size = 8, vjust = 0.7),
          axis.title.y = element_text (size = 8, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 8, vjust = 0.2),
          plot.title = element_text(size = 10, hjust = 0),
          strip.background = element_blank()) +
    ggtitle("B") + 
    scale_x_continuous(breaks = seq(0, 300, 50))
  ,
  
  ggplot(subset(maxvals, Chr != 27), aes(maxcM.male, maxcM.female)) +
    stat_smooth(method = "lm") +
    geom_text(aes(label = Chr), size = 3) +
    annotate("text", x = maxvals$maxcM.male[27], y = maxvals$maxcM.female[27], label = "X", size = 3) +
    coord_cartesian(xlim = c(50, 370)) +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 8, vjust = 0),
          axis.text.y  = element_text (size = 8, hjust = 1.3),
          strip.text.x = element_text (size = 8, vjust = 0.7),
          axis.title.y = element_text (size = 8, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 8, vjust = 0.2),
          plot.title = element_text(size = 10, hjust = 0),
          strip.background = element_blank(),
          legend.position = "top") +
    scale_x_continuous(breaks = seq(0, 400, 50)) +
    scale_y_continuous(breaks = seq(0, 400, 50)) +
    ggtitle("C") + 
    labs(x = "Male Linkage Map Length (cM)", y = "Female Linkage Map Length (cM)")
  ,
  
  cols = 3)
dev.off()

head(maptab)
maptab2 <- select(maptab, SNP.Name, Chr, cMPosition, cMPosition.Female, cMPosition.Male, GenomicPosition)
head(maptab2)
names(maptab2) <- c("SNP.Name", "Chromosome", "cM.Position.Sex.Averaged", "cM.Position.Female",
                    "cM.Position.Male", "Genome.Position.Oarv3.1")

write.table(maptab2, "../Recombination Rate Manuscript/SupplementaryInformation/LinkageMapTable.csv", sep = ",", row.names = F, quote = F)




### PRINT R SQUARED VALUES ON FIGURES

#' Recombination landscape

# crossover probabilities population wide

xoverprob <- read.table("results/2_PopWide_CrossoverPr_bin1000K_Full_g.txt", header = T)
xoverprob <- subset(xoverprob, DistanceToTelo < 61)
xoverprob$MaleToFemale2 <- xoverprob$Male.Pr/xoverprob$Female.Pr
xoverprob$MaleToFemale2 <- xoverprob$Male.Pr-xoverprob$Female.Pr


xoverprob.sex <- droplevels(subset(xoverprob, select = c(Male.Pr, Female.Pr, DistanceToTelo, AdjustedDistance, MaxBin, SNP.Count, BinID, Chr, PC.GC, TelomereClass)))
xoverprob.sex <- melt(xoverprob.sex, id.vars = c("DistanceToTelo", "MaxBin", "SNP.Count", "BinID", "Chr", "AdjustedDistance", "PC.GC", "TelomereClass"))
names(xoverprob.sex)[which(names(xoverprob.sex) %in% c("variable", "value"))] <- c("Sex", "Probability")
xoverprob.sex$Sex <- gsub(".Pr", "", xoverprob.sex$Sex)


sex.interaction.recomb <- glm(Probability ~ Sex * I(DistanceToTelo^3) + Sex * I(DistanceToTelo^2) + Sex*DistanceToTelo + SNP.Count + PC.GC,
                              data = xoverprob.sex,
                              family = gaussian())

summary(sex.interaction.recomb)

sex.interaction.recomb.resids <- glm(Probability ~ SNP.Count + PC.GC, data = xoverprob.sex, family = gaussian())
xoverprob.sex$residuals <- residuals(sex.interaction.recomb.resids)

mf.ratio.recomb <- glm(MaleToFemale ~ I(DistanceToTelo^3) + I(DistanceToTelo^2) + DistanceToTelo + SNP.Count + PC.GC, 
                       data = xoverprob,
                       family = gaussian())

mf.ratio.recomb.resids <- glm(MaleToFemale ~ MaxBin + SNP.Count + PC.GC, data = xoverprob, family = gaussian())
xoverprob$residuals <- residuals(mf.ratio.recomb.resids)

summary(mf.ratio.recomb)

#' polynomial root

model2c <- lm(residuals ~ I(DistanceToTelo^3) + I(DistanceToTelo^2) + DistanceToTelo,data = xoverprob)
tab3 <- summary(model2c)$coefficients
z <- matrix(c(tab3[1,1],tab3[4,1],tab3[3,1],tab3[2,1]), ncol=1) 
polyroot(z)

newpath <- "C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\"

tiff(paste0(newpath, "Fig2.tiff"), width = 5, height = 2.5, units = "in", res = 300) 
multiplot(
 # ggplot(xoverprob.sex, aes(DistanceToTelo, residuals, col = Sex)) +
 ggplot(xoverprob.sex, aes(DistanceToTelo, Probability, col = Sex)) +
    geom_point(alpha = 0.1, size = 1) +
    #stat_smooth(method = "lm", formula = y ~ log(x)) +
    stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x, size = 0.5) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0, 150, 25)) +
    labs(x = "Distance to Nearest Telomere (Mb)", y = "Crossover Probability", col = "") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 8, vjust = 0),
          axis.text.y  = element_text (size = 8, hjust = 1.3),
          strip.text.x = element_text (size = 8, vjust = 0.7),
          axis.title.y = element_text (size = 8, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 8, vjust = 0.2),
          plot.title = element_text(size = 10, hjust = 0),
          strip.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.background = element_blank(),
          legend.position = c(0.8, 0.85))+
    ggtitle("A")
 
  ,
  
  # ggplot(xoverprob, aes(DistanceToTelo, residuals)) +
  ggplot(xoverprob, aes(DistanceToTelo, MaleToFemale)) +
    geom_hline(yintercept = 0) +
    geom_point(alpha = 0.2, size = 1) +
    stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x, size = 0.5) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0, 150, 25)) +
    labs(x = "Distance to Nearest Telomere (Mb)", y = "Male - Female Crossover Probabilities", col = "") +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 8, vjust = 0),
          axis.text.y  = element_text (size = 8, hjust = 1.3),
          strip.text.x = element_text (size = 8, vjust = 0.7),
          axis.title.y = element_text (size = 8, angle = 90, vjust = 1.5),
          axis.title.x = element_text (size = 8, vjust = 0.2),
          plot.title = element_text(size = 10, hjust = 0),
          legend.position = c(0.85, 0.9)) +
    ggtitle("B")
  , cols = 2)
dev.off()


#~~ Output a table of the maximum chromosome values.

maxvals.out <- maxvals[,c(5, 6, 1, 2, 3, 4)]
names(maxvals.out) <- c("Chromosome", "$N_{SNPS}$", "Sex-averaged map (cM)", "Male map (cM)", "Female map (cM)", "Last SNP Position (bp)")

temptable <- xtable(maxvals.out, 
                    label = "ChromosomeStats",
                    align = "ccccccc",
                    caption = paste("Linkage map information for the Soay sheep population using information from",
                                    nrow(recsumm), "sub-pedigrees including", length(unique(recsumm$RRID)),
                                    "unique focal individuals. See Table  \ref{tbl:heritabilitytable} for sex-specific sample sizes."))
print.xtable(temptable, include.rownames = F)



#~~ Output a table of the model results

landscapetab <- rbind(data.frame(summary(sex.interaction.recomb)$coefficients),
                      data.frame(summary(mf.ratio.recomb)$coefficients))
landscapetab$Effect <- row.names(landscapetab)
landscapetab$Model <- ""
landscapetab$Model[grep("Intercept", landscapetab$Effect)] <- c("Crossover Probability", "Male:Female Crossover Probability")

landscapetab <- landscapetab[,c("Model", "Effect", "Estimate", "Std..Error", "t.value", "Pr...t..")]
row.names(landscapetab) <- 1:nrow(landscapetab)

landscapetab <- landscapetab[c(1, 3, 4, 5, 2, 8, 9, 10, 6, 7, 11:16),]
landscapetab$Effect <- c("Female (Intercept)", "Female Dist$^{3}$", "Female Dist$^{2}$", "Female Dist",
                         "Male", "Male Dist$^{3}$", "Male Dist$^{2}$", "Male Dist", "SNP Count", "GC content (%)",
                         "Intercept", "Dist$^{3}$", "Dist$^{2}$", "Dist", "SNP Count", "GC content (%)")

names(landscapetab) <- c("Model", "Effect", "Estimate", "S.E.", "t", "P")

write.table(landscapetab, "../Recombination Rate Manuscript/TabS2_LandscapeModel.txt", sep = "\t", quote = F, row.names = F)

temptable <- xtable(landscapetab, 
                    label = "RecombinationLandscape",
                    align = "ccccccc",
                    caption = paste("Recombination landscape (fix)"))
print.xtable(temptable, include.rownames = F)



# BY GENOTYPE

load("results/2_DataForGenoSpecificRecombLandscapev2.Rdata")

head(genobin.sexgeno)
genobin.fem <- droplevels(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Female", "G.A.Female", "G.G.Female")))
head(genobin.fem)
genobin.male <- droplevels(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Male", "G.A.Male", "G.G.Male")))
head(genobin.male )

geno.model <- rbind(cbind(Sex = "Female",
                                     summary(glm(value ~ variable*I(DistanceToTelo^3) + variable*I(DistanceToTelo^2) + variable*DistanceToTelo + SNP.Count + PC.GC,
                                                 data = genobin.fem))$coefficients),
                               cbind(Sex = "Male",
                                     summary(glm(value ~ variable*I(DistanceToTelo^3) + variable*I(DistanceToTelo^2) + variable*DistanceToTelo + SNP.Count + PC.GC,
                                                 data = genobin.male))$coefficients))
geno.model.effects <- row.names(geno.model)
geno.model <- data.frame(geno.model)

head(geno.model)
for(i in 2:5) geno.model[,i] <- signif(as.numeric(geno.model[,i]), 3)
geno.model <- cbind(Effect = geno.model.effects, geno.model)
names(geno.model) <- c("Model.Term", "Sex", "Effect Size", "S.E.", "t", "P")
geno.model$Model.Term <- gsub("DistanceToTelo", "Dist", geno.model$Model.Term)
geno.model$Model.Term <- gsub("variable", "", geno.model$Model.Term)
geno.model$Model.Term <- gsub("Female", "", geno.model$Model.Term)
geno.model$Model.Term <- gsub("Male", "", geno.model$Model.Term)

write.table(geno.model, "../Recombination Rate Manuscript/GenotypeLandscapeInteraction.txt", row.names = F, quote = F, sep = "\t&\t", eol = "\\\\\n")

pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\CrossoverProbabilityByGenotype.pdf", width = 8, height = 12, useDingbats = F) 
multiplot(
  ggplot(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Female", "G.A.Female", "G.G.Female")),
         aes(DistanceToTelo,  value, col = variable)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x, size = 1) +
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
          plot.title = element_text(size = 18, hjust = 0),
          strip.background = element_blank()) +
    ggtitle("A. Females")
  ,
  ggplot(filter(genobin.sexgeno, DistanceToTelo < 61, variable %in% c("A.A.Male", "G.A.Male", "G.G.Male")),
         aes(DistanceToTelo,  value, col = variable)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x, size = 1) +
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
          plot.title = element_text(size = 18, hjust = 0),
          strip.background = element_blank()) +
  ggtitle("B. Males")
, cols = 1)
dev.off()

#   
#   ggplot(filter(genobin.ratiomelt, DistanceToTelo < 61, value < 50), aes(DistanceToTelo,  value, col = variable)) +
#     geom_point(alpha = 0.2) +
#     stat_smooth(method = "lm", formula = y ~ I(x^3) + I(x^2) + x) +
#     scale_colour_brewer(palette = "Set1") +
#     scale_x_continuous(breaks = seq(0, 150, 25)) +
#     labs(x = "Distance to Nearest Telomere (MB)", y = "Crossover Probability", col = "") +
#     theme_bw() +
#     theme(axis.text.x  = element_text (size = 16, vjust = 0),
#           axis.text.y  = element_text (size = 14, hjust = 1.3),
#           strip.text.x = element_text (size = 16, vjust = 0.7),
#           axis.title.y = element_text (size = 16, angle = 90, vjust = 1.5),
#           axis.title.x = element_text (size = 16, vjust = 0.2),
#           plot.title = element_text(size = 18, hjust = 0),
#           strip.background = element_blank(),
#           legend.background = element_blank())



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 1.1. Distance between double crossovers              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# switchtab.Mid <- read.table("results/2_switchtab.doublexovers.forpaper.txt", header = T, sep = "\t")
# mean.span <- mean(log10(switchtab.Mid$Segment.SpanLength))
# sd.span   <- sd(log10(switchtab.Mid$Segment.SpanLength))
# 
# pdf("../Recombination Rate Manuscript/SupplementaryFigureDoubleXover.pdf", width = 7, height = 6)
# ggplot(switchtab.Mid, aes(log10(Segment.SpanLength))) + 
#   geom_histogram(binwidth = 0.05, col = "black", fill = "grey") +
#   # geom_vline(xintercept = bound95(log10(switchtab.Mid$Segment.SpanLength)), col = "red") +
#   geom_vline(xintercept = mean.span - (2.5*sd.span), linetype = "dashed") +
#   # geom_vline(xintercept = mean.span + (2.5*sd.span)) +
#   scale_fill_brewer(palette = "Set1") +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 12, vjust = 0),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 12, vjust = 0.7),
#         axis.title.y = element_text (size = 12, angle = 90, vjust = 0.2),
#         axis.title.x = element_text (size = 12, vjust = 0.2)) +
#   labs(x = expression(log[10]*" distance (bp) between double-crossovers"))
# dev.off()
# 
# switchtab.Mid$NewColour <- ifelse(switchtab.Mid$Segment.SpanLength < 9.7e6, "black", "darkgrey")
# 
# pdf("../Recombination Rate Manuscript/SupplementaryFigureDoubleXoverv2.pdf", width = 7, height = 6)
# ggplot(switchtab.Mid, aes(Segment.SpanLength, fill = NewColour)) + 
#   geom_histogram(binwidth = 0.5e6) +
#   # geom_vline(xintercept = bound95(log10(switchtab.Mid$Segment.SpanLength)), col = "red") +
#   geom_vline(xintercept = median(switchtab.Mid$Segment.SpanLength), linetype = "dashed") +
#   # geom_vline(xintercept = mean.span + (2.5*sd.span)) +
#   scale_fill_identity() +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 12, vjust = 0),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 12, vjust = 0.7),
#         axis.title.y = element_text (size = 12, angle = 90, vjust = 0.2),
#         axis.title.x = element_text (size = 12, vjust = 0.2)) +
#   labs(x = "Distance (bp) between double-crossovers")
# dev.off()
# 
# 
# mean(switchtab.Mid$Segment.SpanLength)
# median(switchtab.Mid$Segment.SpanLength)

switchtab.Mid <- read.table("results/2_switchtab.doublexovers.forpaperwithSingletons.txt", header = T, sep = "\t")
plottab <-subset(switchtab.Mid, Type == "Mid" & Chr %in% 1:26)

pdf("../Recombination Rate Manuscript/SupplementaryFigureDoubleXoverv2.pdf", width = 7, height = 6)

ggplot(plottab, aes(Segment.SpanLength/1e6, fill = Singleton)) +
  geom_histogram(binwidth = 2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Span length of double crossover (Mb)", fill = "Single SNP Crossover?") +
  geom_vline(xintercept = c(9.7)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 400, 50)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "top")

dev.off()

table(plottab$Singleton, plottab$Segment.SpanLength < 9.7e6)


switchtab.Mid$NewColour <- ifelse(switchtab.Mid$Segment.SpanLength < 9.7e6, "black", "darkgrey")

pdf("../Recombination Rate Manuscript/SupplementaryFigureDoubleXoverv2.pdf", width = 7, height = 6)
ggplot(switchtab.Mid, aes(Segment.SpanLength, fill = NewColour)) + 
  geom_histogram(binwidth = 0.5e6) +
  # geom_vline(xintercept = bound95(log10(switchtab.Mid$Segment.SpanLength)), col = "red") +
  geom_vline(xintercept = median(switchtab.Mid$Segment.SpanLength), linetype = "dashed") +
  # geom_vline(xintercept = mean.span + (2.5*sd.span)) +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, vjust = 0),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12, vjust = 0.7),
        axis.title.y = element_text (size = 12, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 12, vjust = 0.2)) +
  labs(x = "Distance (bp) between double-crossovers")
dev.off()


mean(switchtab.Mid$Segment.SpanLength)
median(switchtab.Mid$Segment.SpanLength)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 2. Heritability of recombination rate                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# heritability results

load("results/2_GRMheritability_results_for_paper.txt")




heritability.results$grm.fixef
heritability.results$grm.ranef
heritability.results$bivar.fixef
heritability.results$bivar.ranef
heritability.results$bivar.sig.diff.from.zero
heritability.results$bivar.sig.diff.from.one
heritability.results$bivar.sig.diff.Va
heritability.results$bivar.sig.diff.Vr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 3. Genetic architecture of autosomal recomb rate     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Chromosome partitioning

chr.partitions <- rbind(read.table("results/2_Chromosome_Partition_Logliks.txt", sep = "\t", header = T, stringsAsFactors = F),
                        read.table("results/2_Chromosome_Partition_LogliksCistrans.txt", sep = "\t", header = T, stringsAsFactors = F))

head(chr.partitions)

chr.partitions$Pval <- 1 - pchisq(2*(chr.partitions$LogL - chr.partitions$LogL.wochr), df = 1)
chr.partitions$sig <- ifelse(chr.partitions$Pval < 0.01, "**", ifelse(chr.partitions$Pval < 0.05, "*", ""))


chreffect <- subset(chr.partitions, VarComp == "giv(RRID2).giv")

chreffect$Pval <- 1 - pchisq(2*(chreffect$LogL - chreffect$LogL.wochr), df = 1)
chreffect$sig <- ifelse(chreffect$Pval < 0.01, "**", ifelse(chreffect$Pval < 0.05, "*", ""))


chreffect$SEX[which(chreffect$SEX == "All")] <- "Both Sexes"
chreffect$Label <- NA

comparison.lms <- data.frame(SEX = rep(sort(unique(chreffect$SEX)), times = 2),
                             Analysis = rep(sort(unique(chreffect$Analysis)), each = 3),
                             Prefix = rep(c("a", "b", "c"), times = 2),
                             Label = NA)

comparison.lms$AdjR2 <- NA
comparison.lms$P <- NA

for(i in 1:nrow(comparison.lms)){
  fitvec <- which(chreffect$SEX == comparison.lms$SEX[i] & chreffect$Analysis == comparison.lms$Analysis[i])
  fit <- summary(lm(Effect ~ maxGenome, data = chreffect[fitvec,]))
  comparison.lms$AdjR2[i] <- fit$adj.r.squared
  comparison.lms$P[i] <- fit$coefficients[2,4]
  chreffect$Label[fitvec] <- paste0(comparison.lms$Prefix[i], ". ", comparison.lms$SEX[i])
  comparison.lms$Label[i] <- paste0(comparison.lms$Prefix[i], ". ", comparison.lms$SEX[i])
}

comparison.lms$AdjR2 <- round(comparison.lms$AdjR2, digits = 3)
comparison.lms$P <- round(comparison.lms$P, digits = 3)
# 
# pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\ChromosomePartitioning.pdf", width = 15, height = 5) 
# 
# ggplot(subset(chreffect, Analysis == "trans only"), aes(maxGenome, Effect)) +
#   stat_smooth(method = "lm") +
#   geom_hline(yintercept = 0, colour = "darkgrey") +
#   geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0) +
#   geom_point(alpha = 1, size = 8, col = "black") +
#   geom_text(data=subset(comparison.lms, Analysis == "trans only"), aes(x=220000000, y=0.15, label=paste("Adj. R2 =", AdjR2, "\nP = ", P))) +
#   #geom_point(alpha = 1, size = 8, col = "white") +
#   geom_text(aes(label = Chr), size = 4, col = "white", fontface = "bold") +
#   geom_text(aes(label = sig), size = 8, hjust=-0.5, vjust=0) +
#   facet_wrap(~Label) +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 12),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 16, vjust = 0.7, hjust = 0),
#         axis.title.y = element_text (size = 14, vjust = 0.9),
#         axis.title.x = element_text (size = 14, vjust = 0.2),
#         strip.background = element_blank()) +
#   scale_x_continuous(breaks = seq(0, max(chreffect$maxGenome), 50000000), labels = seq(0, max(chreffect$maxGenome), 50000000)/1e6) +
#   labs(x = "Chromosome Length (MB)", y = "Proportion of VP explained")
# 
# dev.off()

pdf("C:\\Users\\Susan Johnston\\Desktop\\Recombination Rate Manuscript\\ChromosomePartitioningCisTrans.pdf", width = 15, height = 5, useDingbats = F) 

ggplot(subset(chreffect, Analysis == "cis&trans"), aes(maxGenome, Effect)) +
  stat_smooth(method = "lm") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0) +
  geom_point(alpha = 1, size = 8, col = "black") +
  #geom_point(alpha = 1, size = 8, col = "white") +
  geom_text(aes(label = Chr), size = 4, col = "white", fontface = "bold") +
  geom_text(aes(label = sig), size = 8, hjust=-0.5, vjust=0) +
  geom_text(data=subset(comparison.lms, Analysis == "cis&trans"), aes(x=220000000, y=0.16, label=paste("Adj. R2 =", AdjR2, "\nP = ", P))) +
  facet_wrap(~SEX) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank())+
  scale_x_continuous(breaks = seq(0, max(chreffect$maxGenome), 50000000), labels = seq(0, max(chreffect$maxGenome), 50000000)/1e6) +
  labs(x = "Chromosome Length (MB)", y = "Proportion of VP explained")

dev.off()


chrpaper <- filter(chr.partitions, Analysis == "cis&trans", VarComp == "giv(RRID2).giv")
chrpaper <- select(chrpaper, component, std.error, Effect, SE, Chr, SEX, maxGenome, nsnps, Pval)

castab <- cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "component")
names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.VA", "Female.VA", "Male.VA")
castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "std.error"))
names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.VAse", "Female.VAse", "Male.VAse")
castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "Effect"))
names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.h2", "Female.h2", "Male.h2")
castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "SE"))
names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.h2se", "Female.h2se", "Male.h2se")
castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "Pval"))
names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.Pval", "Female.Pval", "Male.Pval")

head(castab)
castab <- arrange(castab, Chr)

for(i in 4:ncol(castab)) castab[,i] <- round(castab[,i], 3)

castab <- castab[,c(3, 1, 2, 4, 7, 10, 13, 16, 5, 8, 11, 14, 17, 6, 9, 12, 15, 18)]
castab$maxGenome <- castab$maxGenome/1e6
castab$maxGenome <- round(castab$maxGenome, 1)

for(i in 4:ncol(castab)) castab[,i] <- round(castab[,i], 3)
test <- colSums(castab)

for(i in c(5, 7, 10, 12, 15, 17)){
  castab[,i] <- paste0("(", castab[,i], ")")
  castab[,i-1] <- paste(castab[,i-1], castab[,i])
}
castab <- rbind(castab, test)

castab <- castab[,which(!1:ncol(castab) %in% c(5, 7, 10, 12, 15, 17))]

write.table(castab, "../Recombination Rate Manuscript/ChromosomePartitioningCisTrans.txt", sep = " & ", quote = F, row.names = F, eol = " \\\\\n")

###
# chrpaper <- filter(chr.partitions, Analysis == "trans only", VarComp == "giv(RRID2).giv")
# chrpaper <- select(chrpaper, component, std.error, Effect, SE, Chr, SEX, maxGenome, nsnps, Pval)
# 
# castab <- cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "component")
# names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.VA", "Female.VA", "Male.VA")
# castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "std.error"))
# names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.VAse", "Female.VAse", "Male.VAse")
# castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "Effect"))
# names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.h2", "Female.h2", "Male.h2")
# castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "SE"))
# names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.h2se", "Female.h2se", "Male.h2se")
# castab <- join(castab, cast(chrpaper, maxGenome + nsnps + Chr ~ SEX, value = "Pval"))
# names(castab)[which(names(castab) %in% c("All", "Female", "Male"))] <- c("All.Pval", "Female.Pval", "Male.Pval")
# 
# head(castab)
# castab <- arrange(castab, Chr)
# 
# for(i in 4:ncol(castab)) castab[,i] <- round(castab[,i], 3)
# 
# castab <- castab[,c(3, 1, 2, 4, 7, 10, 13, 16, 5, 8, 11, 14, 17, 6, 9, 12, 15, 18)]
# castab$maxGenome <- castab$maxGenome/1e6
# castab$maxGenome <- round(castab$maxGenome, 1)
# 
# for(i in 4:ncol(castab)) castab[,i] <- round(castab[,i], 3)
# test <- colSums(castab)
# 
# for(i in c(5, 7, 10, 12, 15, 17)){
#   castab[,i] <- paste0("(", castab[,i], ")")
#   castab[,i-1] <- paste(castab[,i-1], castab[,i])
# }
# castab <- rbind(castab, test)
# 
# castab <- castab[,which(!1:ncol(castab) %in% c(5, 7, 10, 12, 15, 17))]
# 
# write.table(castab, "../Recombination Rate Manuscript/ChromosomePartitioning.txt", sep = " & ", quote = F, row.names = F, eol = " \\\\\n")




#' regional heritability

analysis.vec <- c("Run_a.ps_150SNP_ChrChrall_trans",
                  "Run_b.ps_50SNP_ChrChrall_trans",
                  "Run_c.ps_20SNP_ChrChrall_trans",
                  "Run_e.ps_10SNP_Chr6e_trans",
                  "Run_h.ps_10SNP_Chr6end_trans",
                  "Run_h.ps_5SNP_Chr6end_trans")




snplistvec.150 <- read.table(paste0("results/2_Regional_Heritability_Run_a.ps_150SNP_Chrall_cistrans.txt"), header = T)
snplistvec.150$Chr <- factor(snplistvec.150$Chr)
snplistvec.150$Colour <- factor(snplistvec.150$Colour)
snplistvec.150$Label <- recoderFunc(snplistvec.150$Sex,
                                    c("Both", "Female", "Male"),
                                    c("a. Both Sexes", "b. Female", "c. Male"))

# ggplot(snplistvec.150, aes(MedianCumuPos, Reg.h2, col = Colour, group = Chr)) + 
#   geom_point(size = 2, aes(shape = Error)) +
#   geom_errorbar(aes(ymin =  Reg.h2-Reg.h2.SE, ymax = Reg.h2+Reg.h2.SE), col = "grey", alpha = 0.3) +
#   theme_bw() +
#   geom_line() +
#   scale_color_brewer(palette = "Set1") +
#   theme(axis.text.x  = element_text (size = 12),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
#         axis.title.y = element_text (size = 14, vjust = 0.9),
#         axis.title.x = element_text (size = 14, vjust = 0.2),
#         strip.background = element_blank()) +
#   theme(legend.position = "none") +
#   scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
#                      labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
#   labs(x = "Chromosome", y = "-log10(P)") +
#   facet_wrap(~Label, ncol = 1)

pdf("../Recombination Rate Manuscript/150snpwindowsRh2.pdf", width = 8, height = 8)
ggplot(snplistvec.150, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2, aes(shape = Error)) +
  theme_bw() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.150)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()

arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6)), Reg.h2.constraint == "Positive"), Sex, MedianCumuPos)
arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6)), Chr == 1), Sex)


snplistvec.50 <- read.table(paste0("results/2_Regional_Heritability_Run_b.ps_50SNP_Chrall_cistrans.txt"), header = T)
snplistvec.50$Chr <- factor(snplistvec.50$Chr)
snplistvec.50$Colour <- factor(snplistvec.50$Colour)
snplistvec.50$Label <- recoderFunc(snplistvec.50$Sex,
                                    c("Both", "Female", "Male"),
                                    c("a. Both Sexes", "b. Female", "c. Male"))


pdf("../Recombination Rate Manuscript/50snpwindowsRh2.pdf", width = 8, height = 8)
ggplot(snplistvec.50, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2, aes(shape = Error)) +
  theme_bw() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.50)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()


snplistvec.20 <- read.table(paste0("results/2_Regional_Heritability_Run_c.ps_20SNP_Chrall_cistrans.txt"), header = T)
snplistvec.20$Chr <- factor(snplistvec.20$Chr)
snplistvec.20$Colour <- factor(snplistvec.20$Colour)
snplistvec.20$Label <- recoderFunc(snplistvec.20$Sex,
                                   c("Both", "Female", "Male"),
                                   c("A. All Sheep", "B. Females", "C. Males"))
newpath <- "C:/Users/Susan Johnston/Desktop/Drive/Journal Submissions/12. Johnston RNF212 Recombination Paper/PLoS Genetics/First Submission/"

tiff(paste0(newpath, "Fig3.tiff"), width = 5, height = 5, units = "in", res = 300) 
#pdf("../Recombination Rate Manuscript/20snpwindowsRh2.pdf", width = 8, height = 8, useDingbats = F)
ggplot(snplistvec.20, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 1) +
  #geom_point(size = 2, aes(shape = Error)) +
  #geom_line(data = snplistvec.150, aes(MedianCumuPos, -log10(P)), alpha = 0.5, colour = "grey") + 
  theme_bw() +
  geom_line(size = 0.4) +
  scale_color_manual(values = c("grey20", "grey55")) +
  theme(axis.text.x  = element_text (size = 8),
        axis.text.y  = element_text (size = 8),
        strip.text.x = element_text (size = 8, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 8, vjust = 0.9),
        axis.title.x = element_text (size = 8, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.20)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()



ggplot(snplistvec.20, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2) +
  #geom_point(size = 2, aes(shape = Error)) +
  #geom_line(data = snplistvec.150, aes(MedianCumuPos, -log10(P)), alpha = 0.5, colour = "grey") + 
  theme_bw() +
  geom_line(size = 0.4) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 20),
        axis.text.y  = element_text (size = 20),
        strip.text.x = element_text (size = 20, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 20, vjust = 0.9),
        axis.title.x = element_text (size = 20, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.20)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)







head(snplistvec.20)

res150 <- arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6))), Sex)
res50 <- arrange(filter(snplistvec.50, P < (0.05/(nrow(snplistvec.50)/6))), Sex)
res20 <- arrange(filter(snplistvec.20, P < (0.05/(nrow(snplistvec.20)/6))), Sex)


filter(snplistvec.20,  V1 == "chr6.c.ps.1827.20")


regh2.results <- rbind(res150, res50, res20)
regh2.results <- filter(regh2.results, snpcount %in% c(20, 50, 150))
regh2.results <- select(regh2.results, V1, Sex, Chr, FirstPos, LastPos, Reg.h2, Reg.h2.SE, Genome.Rest.h2, Genome.Rest.h2.SE, P, Error)
filter(regh2.results, Reg.h2 > 0.01)


fullregh2 <- rbind(cbind(Window = 150, snplistvec.150),
                   cbind(Window = 50, snplistvec.50),
                   cbind(Window = 20, snplistvec.20))

head(fullregh2)
fullregh2 <- select(fullregh2,Window, Sex,  FirstPos, LastPos, Chr,  
                    Genome.Rest.h2, Genome.Rest.h2.SE, Reg.h2, Reg.h2.SE, 
                    Error, Reg.h2.constraint, Genome.GRM.LL, Genome.Reg.GRM.LL, Chi.Value, P)

names(fullregh2) <- c("Window.Size", "Sex", "First.SNP.Pos", "Last.SNP.Pos", "Chromosome",
                      "RestofGenome.GRM.h2", "RestofGenome.GRM.h2.SE", 
                      "Window.GRM.h2", "Window.GRM.h2.SE", "ASReml.Error",
                      "Window.GRM.h2.constraint", "RestofGenome.GRM.LL", 
                      "RestofGenome.and.Window.GRM.LL", "Chi.Value", "P")



write.table(fullregh2, "../Recombination Rate Manuscript/SupplementaryInformation/RegionalHeritabilityResults.csv", row.names = F, sep = ",", quote = F) 
  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 4. Genome wide association studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# GWAS results cis&trans

gwas.cistrans  <- read.table("results/2_GWAS_Results_cistrans.txt", header = T, stringsAsFactors = F)
gwas.transonly <- read.table("results/2_GWAS_Results_trans.txt", header = T, stringsAsFactors = F)

gwas.cistrans <- join(gwas.cistrans  , select(maptab, SNP.Name, Cumu2))
gwas.transonly <- join(gwas.transonly, select(maptab, SNP.Name, Cumu2))

gwas.cistrans$Cumu2[which(is.na(gwas.cistrans$Cumu2))] <- 0
gwas.transonly$Cumu2[which(is.na(gwas.transonly$Cumu2))] <- 0

gwas.cistrans$Chromosome[which(gwas.cistrans$Chromosome == "X")] <- "27"
gwas.transonly$Chromosome[which(gwas.transonly$Chromosome == "X")] <- "27"

gwas.cistrans$Chromosome <- as.numeric(as.character(gwas.cistrans$Chromosome))
gwas.transonly$Chromosome <- as.numeric(as.character(gwas.transonly$Chromosome))

gwas.cistrans$Analysis[which(gwas.cistrans$Analysis == "BothSexes")] <- "Both Sexes"
gwas.transonly$Analysis[which(gwas.transonly$Analysis == "all")] <- "Both Sexes"
gwas.transonly$Analysis[which(gwas.transonly$Analysis == "female")] <- "Female"
gwas.transonly$Analysis[which(gwas.transonly$Analysis == "male")] <- "Male"


gwas.cistrans$Colour <- "1"
gwas.cistrans$Colour[which(gwas.cistrans$Chr %in% seq(0, 26, 2))] <- "2"
gwas.transonly$Colour <- "1"
gwas.transonly$Colour[which(gwas.transonly$Chr %in% seq(0, 26, 2))] <- "2"


# Thresholds and plot info

lambda <- tapply(gwas.transonly$chi2.1df, gwas.transonly$Analysis, median)/median(rchisq(nrow(gwas.transonly)/3, 1))
Pthresh <- 0.05/22273.61


chrinfo <- NULL

for(i in unique(gwas.cistrans$Chromosome)){
  
  temp1 <- subset(gwas.cistrans, Chromosome == i)
  temp2 <- data.frame(Chromosome = i,
                      Start = temp1[1,"Cumu2"],
                      Stop = temp1[nrow(temp1),"Cumu2"])
  chrinfo <- rbind(chrinfo, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$Chromosome2 <- c(0:10,"", 12, "", 14, "", 16, "", 18, "", "", 21, "",  "", 24, "", "", "X")


# Plot


gwas.cistrans$Label <- recoderFunc(gwas.cistrans$Analysis,
                                   c("Both Sexes", "Females", "Males"),
                                   c("A. Both Sexes", "B. Females", "C. Males"))

gwas.transonly$Label <- recoderFunc(gwas.transonly$Analysis,
                                   c("Both Sexes", "Females", "Males"),
                                   c("A. Both Sexes", "B. Females", "C. Males"))


gwas <- ggplot(gwas.cistrans, aes(Cumu2, -log10(Pc1df), col = factor(Colour))) +
  geom_point(size = 2.5,alpha = 0.7) +
  geom_hline(yintercept = -log10(Pthresh),linetype = 2, alpha = 0.6, size = 1) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = chrinfo$Mid,
                     labels = c(0:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26, "X")) +
  scale_y_continuous(breaks = seq(0, 14, 1)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
gwas

pdf("../Recombination Rate Manuscript/GWAScistrans.pdf", width = 8, height = 8, useDingbats = F)
gwas
dev.off()


gwa_res.null<-NULL

for(i in sort(unique(gwas.cistrans$Label))){
  subdata <- subset(gwas.cistrans, Label == i)
  x <- data.frame(obs = sort(-log10(subdata$Pc1df)),
                  exp = sort(-log10(seq(1/nrow(subdata),1,1/nrow(subdata)))),
                  Label = i)
  gwa_res.null <- rbind(gwa_res.null, x)
  rm(subdata, x)
}
gc()

tapply(gwas.cistrans$chi2.1df, gwas.cistrans$Label, median)

tapply(gwas.cistrans$chi2.1df, gwas.cistrans$Label, median)/median(rchisq(nrow(gwas.cistrans)/3, 1))

ppplot <- ggplot(gwa_res.null, aes(exp, obs)) +
  geom_point(size = 2.5,alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1,linetype = 2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0, colour = "white"),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 14, 1)) +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
  facet_wrap(~Label, ncol = 1)


pdf("../Recombination Rate Manuscript/GWAScistrans.pdf", width = 10, height = 8, useDingbats = F)
multiplot(gwas, ppplot, cols = 2, layout = matrix(c(1, 1, 2), nrow = 1))
dev.off()


gwas.cistrans[which(gwas.cistrans$Pc1df < Pthresh),]


gwastab.supp <- select(gwas.cistrans, -Cumu2, -Colour, -Label)
gwastab.supp <- select(gwastab.supp, Analysis, SNP.Name, Chromosome, Position, A1, A2, N, CallRate, MAF, effAB, effBB, Pc1df)

head(gwastab.supp)
names(gwastab.supp) <- c("Sex", "SNP.Name", "Chromosome", "Position.bp", "Allele1", "Allele2", "N", "Call.Rate", "Minor.Allele.Freq",
                         "EffectAB", "EffectBB", "P.Corrected")

write.table(gwastab.supp, "../Recombination Rate Manuscript/SupplementaryInformation/GWASCisTransResults.csv", row.names = F, sep = ",", quote = F) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Determine the proportion of additive genetic variance explained     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ Get the results of an animal model with s74824.1 genotype fitted as a fixed effect

load("results/2_rnf_cistransmodels.Rdata", verbose = T)

temp <- row.names(rnf.cistrans.fixef)
effect.snp <- data.frame(rnf.cistrans.fixef)
for(i in 1:3) effect.snp[,i] <- as.numeric(as.character(effect.snp[,i]))
effect.snp$Effect <- temp
effect.snp <- dplyr::filter(effect.snp, Analysis != "Both Sexes")
effect.snp2 <- effect.snp[grep("Intercept", effect.snp$Effect),]

effect.snp <- effect.snp[grep("s748", effect.snp$Effect),]
effect.snp$Effect <- rep(c("A/A", "G/A", "G/G"), times = 2)
effect.snp

effect.snp$solution2 <- effect.snp$solution
effect.snp$solution2[which(effect.snp$Effect == "A/A")] <- effect.snp2$solution
effect.snp$solution2se <- effect.snp$std.error
effect.snp$solution2se[which(effect.snp$Effect == "A/A")] <- effect.snp2$std.error

effect.snp$solution3 <- NA
effect.snp$solution3[1:3] <- effect.snp$solution[1:3] + effect.snp2$solution[1]
effect.snp$solution3[4:6] <- effect.snp$solution[4:6] + effect.snp2$solution[2]


sigsnptab <- data.frame(SigGeno = as.character.gwaa.data(sheepabel[,c("s74824.1", "s10844.1", "s14138.1")]))
sigsnptab$RRID <- row.names(sigsnptab)
recsumm <- join(recsumm, sigsnptab)
rectab <- join(rectab, sigsnptab)


data.frame(table(recsumm$RRID.SEX, recsumm$SigGeno.s74824.1))

table(recsumm$SigGeno.s74824.1, recsumm$SigGeno.s14138.1, useNA = "always")
table(recsumm$SigGeno.s74824.1, recsumm$SigGeno.s10844.1, useNA = "always")
table(recsumm$SigGeno.s14138.1, recsumm$SigGeno.s10844.1, useNA = "always")


ggplot(recsumm, aes(SigGeno.s74824.1, TotalRecombCount)) + 
  geom_boxplot(notch = T) + 
  facet_wrap(~RRID.SEX)+
  labs(x = "s74824.1 Genotype", y = "Total Recombination Count") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank())



effect.snp$Label <- rep(c("B. Males", "A. Females"), each = 3)

temp <- data.frame(table(recsumm$RRID.SEX, recsumm$SigGeno.s74824.1))
names(temp)[3] <- "Observations"

temp2 <- unique(subset(recsumm, select = c(RRID, RRID.SEX, SigGeno.s74824.1)))

temp3 <- data.frame(table(temp2$RRID.SEX, temp2$SigGeno.s74824.1))
names(temp3)[3] <- "UniqueIDs"

temp <- join(temp, temp3)
names(temp)[1:2] <- c("Analysis", "Effect")
temp$Analysis <- paste0(temp$Analysis, "s")

effect.snp <- join(effect.snp, temp)
rm(temp, temp2, temp3)


effect.snp.m <- subset(effect.snp, Analysis == "Males")
effect.snp.f <- subset(effect.snp, Analysis == "Females")

pdf("../Recombination Rate Manuscript/snpeffectsizes.pdf", width = 6, height = 4, useDingbats = F)
#postscript("../Recombination Rate Manuscript/snpeffectsizes.ps", width = 600, height = 400)

multiplot(
ggplot(effect.snp.f, aes(Effect, solution3)) + 
  geom_rect(xmin = 0, xmax = 4, ymin = effect.snp.f$solution3[1] + 0.8, ymax = effect.snp.f$solution3[1] + 2, fill = "white", col = "black") +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymax = solution3 + solution2se, ymin = solution3 - solution2se), width = 0) +
  facet_wrap(~Label) +
  labs(x = "s74824.1 Genotype", y = "Effect size") +
  theme_bw() +
  coord_cartesian(ylim = c(effect.snp.f$solution3[1] - 5, effect.snp.f$solution3[1] +2)) +
  geom_text(aes(y = effect.snp.f$solution3[1] + 1.7, label = Observations), size = 4) +
  geom_text(aes(y = effect.snp.f$solution3[1] + 1.2, label = paste0("(", UniqueIDs, ")")), size = 4) +
  scale_y_continuous(breaks = seq(0, 29, 1)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank())
,
ggplot(effect.snp.m, aes(Effect, solution3)) + 
  geom_rect(xmin = 0, xmax = 4, ymin = effect.snp.m$solution3[1] + 0.8, ymax = effect.snp.m$solution3[1] + 2, fill = "white", col = "black") +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymax = solution3 + solution2se, ymin = solution3 - solution2se), width = 0) +
  facet_wrap(~Label) +
  labs(x = "s74824.1 Genotype", y = "Effect size") +
  theme_bw() +
  coord_cartesian(ylim = c(effect.snp.m$solution3[1] - 5, effect.snp.m$solution3[1] +2)) +
  geom_text(aes(y = effect.snp.m$solution3[1] + 1.7, label = Observations), size = 4) +
  geom_text(aes(y = effect.snp.m$solution3[1] + 1.2, label = paste0("(", UniqueIDs, ")")), size = 4) +
  scale_y_continuous(breaks = seq(0, 35, 1)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank())
,cols = 2)

dev.off()

beepr::beep()




ggplot(effect.snp, aes(Effect, solution)) + 
  geom_rect(xmin = 0, xmax = 4, ymin = 0.01, ymax = 1.5, fill = "white") +
  geom_text(aes(y = 1, label = Observations), size = 4) +
  geom_text(aes(y = 0.5, label = paste0("(", UniqueIDs, ")")), size = 4) +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymax = solution + std.error, ymin = solution - std.error), width = 0) +
  facet_wrap(~Label, scales = "free_y") +
  labs(x = "s74824.1 Genotype", y = "Effect size (relative to intercept at 0)") +
  theme_bw() +
  scale_y_continuous(breaks = seq(-4, 0, 1)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 1.4),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank())



effect.snp.latex <- select(effect.snp,Analysis,Effect,Observations,UniqueIDs,solution2,solution2se,z.ratio)
names(effect.snp.latex) <- c("Sex", "s74824.1 Genotype", "Nobs",
                             "Nids", "Effect Size", "S.E.", "Z Ratio")
write.table(effect.snp.latex, "../Recombination Rate Manuscript/SNPeffectTable.txt",
            sep = " & ", eol = "\\\\\n", row.names = F, quote = F)


#~~ Load library
# {
# library(msm)
# 
# #~~ Specify values
# 
# # NB. We adopt the convention that the A allele is the allele which increases
# # the trait value. In Cyp study, genotype 3 appears to increase the trait value.
# # Therefore the levels 1, 2 and 3 for CYP should correspond to genotypes BB, AB
# # and AA, respectively.
# 
# p <- 0.4264862                # Frequency of A allele
# q <- 0.5735138               # Frequency of B allele
# 
# # Specify fixed effect sizes as estimated from the animal model (.asr file)
# 
# effectAA <- 0             # Effect Size of genotype AA
# effectAB <- -1.58932907884297       # Effect Size of genotype AB
# effectBB <- -3.27439771041814      # Effect Size of genotype BB
# 
# # Specify the additive genetic variance from same animal model (.asr file)
# 
# Va <-  4.082639               # Additive Genetic Variance
# 
# # Specify the variance covariance matrix of the fixed effects. (.vrb file)
# # e.g. for the email example,
# # CYP2     0.882053E+01
# # CYP1     0.563167E+01   0.161131E+02
# 
# VarAB <- 0.1860764
# VarBB <- 0.2435467 
# CovAB.BB <- 0.1594847
# 
# 
# a <- (effectAA - effectBB)/2     # Calculate a
# d <-  a + effectAB               # Calculate d
# 
# Vq <- 2*p*q*(a + d*(q - p))^2        # Calculate Vq
# VarExplained <- Vq/(Vq + Va)         # Calculate additive genetic variance explained by QTL
# 
# 
# vcovMatrix <- matrix(data=c(VarAB, CovAB.BB, CovAB.BB, VarBB), nrow=2)   # VarAB, CovABBB, CovABBB, VarBB
# 
# 
# x1 <- vcovMatrix[1,1]   #  variance of the het effect
# x2 <- vcovMatrix[2,2]   #  variance of the hom effect
# 
# beta <- c(effectAB, effectBB)   # create vactor beta for deltamethod
# 
# X <- 2*p*q
# Y <- q^2
# 
# Vq.se <- deltamethod(~X*(-x2/2 + (-x2/2 + x1)*Y)^2, beta, vcovMatrix)  # standard error
# 
# results <- list(Vq, Va, VarExplained, Vq.se)
# names(results) <- c("Vsnp", "Va", "Va Explained by Vsnp", "Vsnp SE")
# 
# results
# 
# # system("cmd", input = "cp gcta/Run_c.ps_20SNP_Chrall_cistrans/chr6.c.ps.1827.20snplist.txt data/chr6.c.ps.1827.20snplist.txt")
# 
# # plink --file 20150129merged1_66nodups.QC2 --from s70829.1 --to s74824.1 --sheep --recode --out chr6segmentmach
# # plink --file 20140214_SheepHD_QC1_Polym --from s70829.1 --to oar3_OAR6_117013262 --sheep --recode --out chr6HDsegmentmach
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Gene information                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
#' Add gene information

# Gene positions

genes <- read.table("data/Ovis_aries.Oar_v3.1.79.genes.txt", header = T, sep = "\t", stringsAsFactors = F)
exons <- read.table("data/Ovis_aries.Oar_v3.1.79.exons.txt", header = T, sep = "\t", stringsAsFactors = F)
exons <- filter(exons, Chromosome == 6, Start > 1e8)
gc()

genes <- rbind(genes, c(6, "unknown", "gene", 116430959, 116442905, "-", "unknown", "RNF212", "protein_coding"))       
genes$Start <- as.numeric(genes$Start)
genes$Stop  <- as.numeric(genes$Stop)

genes$MidPos <- apply(genes[,c("Start", "Stop")], 1, mean)
genes <- arrange(genes, Chromosome, Start)
genes$Order <- 1:nrow(genes)
genes$GeneLength <- genes$Stop - genes$Start

#~~ Chromosome 6
genes2 <- filter(genes, Chromosome == 6, Start > 1e8)
genes2 <- filter(genes2, !gene_id %in% c("ENSOARG00000011847", "ENSOARG00000015218"), gene_biotype == "protein_coding")
genes2$newCode <- rep(1:8, length.out = nrow(genes2))
genes2$newCode[which(genes2$gene_name == "LRPAP1")] <- 7
genemelt <- melt(select(genes2, Chromosome, Start, Stop, gene_name, Order, MidPos, GeneLength, newCode),
                 id.vars = c("Chromosome", "gene_name", "Order", "MidPos", "GeneLength", "newCode"))

nrow(genes2[which(genes2$Stop > 115207925),])

genes3 <- filter(genes, Chromosome == 7, Stop > 20266922, Start < 21355824, gene_biotype == "protein_coding")



scangraph <- filter(snplistvec.20, Chr == 6, Sex == "Female")
scangraph <- scangraph[-nrow(scangraph),]
scangraph2 <- melt(select(scangraph, FirstPos, LastPos, MedianPos, Chr, P, Colour, Error, Order),
                   id.vars = c("MedianPos", "Chr", "P", "Colour", "Error", "Order"))

gwasgraph <- filter(gwas.cistrans, Chromosome == 6, Analysis == "Females")


ggplot() +
  geom_rect(aes(xmax = 118e8, xmin = 0,  ymax = 23, ymin = 14), fill = "white", col = "black") +
  geom_line(data = scangraph2, aes(value,    -log10(P), group = Order), size = 3,  alpha = 0.3) +
  geom_line(data = genemelt  , aes(value,    -log10(min(scangraph$P, na.rm = T)) + newCode + 1, group = Order), size = 2) +
  geom_text(data = genes2    , aes(MidPos,   -log10(min(scangraph$P, na.rm = T)) + newCode + 1.4, label = gene_name), size = 3) +
  geom_point(data = gwasgraph, aes(Position, -log10(Pc1df)), size = 4, alpha = 0.5) +
  coord_cartesian(xlim = c(1.15e8-50000, 1.175e8+50000), ylim = c(-0.5, 23)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.9, hjust = 0.3),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.20)/6))) + 
  geom_hline(yintercept = -log10(Pthresh), linetype = "dashed") + 
  geom_hline(yintercept = -log10(Pthresh), linetype = "dashed") + 
  labs(x = "Base Pair Position (Mb)", y = "-log10(P)") +
 scale_x_continuous(breaks = seq(0, 120e6, 2.5e6), labels = seq(0, 120, 2.5)) +
 scale_y_continuous(breaks = seq(0, 12, 2), labels = seq(0, 12, 2))

#~~ Mach analysis

chr6results <- read.table("mach/ImputationDataForPlotting.txt", header = T)

ggplot() + 
  geom_rect(aes(xmax = 118e8, xmin = 0,  ymax = 23, ymin = 14), fill = "white", col = "black") +
  geom_text(data = genes2, aes(x = MidPos, y = 9, label = gene_name), angle = 270, size = 4) +
  geom_line(data = genemelt, aes(value, 8, group = Order), size = 6) +
  geom_point(data = chr6results, aes(Position, -log10(Pc1df), colour = Array, shape = Array), size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(Pthresh),linetype = 2, alpha = 0.6, size = 1) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 118000000, 200000), labels = seq(0, 118, 0.2)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x = "Position (MB)", y = "-log10(P)") +
  coord_cartesian(ylim = c(-0.5, 10), xlim = c(115700000, max(chr6results$Position) + 0.5e5)) +
  scale_y_continuous(breaks = seq(0, 12, 2), labels = seq(0, 12, 2))

# 
# ggplot() +
#   geom_rect(aes(xmax = 118e8, xmin = 0,  ymax = 23, ymin = 14), fill = "white", col = "black") +
#   geom_line(data = scangraph2, aes(value,    -log10(P), group = Order), size = 3,  alpha = 0.3) +
#   geom_line(data = genemelt  , aes(value,    -log10(min(scangraph$P, na.rm = T)) + newCode + 1, group = Order), size = 2) +
#   geom_text(data = genes2    , aes(MidPos,   -log10(min(scangraph$P, na.rm = T)) + newCode + 1.4, label = gene_name), size = 3) +
#   geom_point(data = gwasgraph, aes(Position, -log10(Pc1df)), size = 4, alpha = 0.5) +
#   coord_cartesian(xlim = c(1.15e8-50000, 1.175e8+50000), ylim = c(-0.5, 23)) +
#   theme_bw() +
#   theme(axis.text.x  = element_text (size = 16, vjust = 0),
#         axis.text.y  = element_text (size = 14),
#         strip.text.x = element_text (size = 16, vjust = 0.7),
#         axis.title.y = element_text (size = 16, angle = 90, vjust = 0.9, hjust = 0.3),
#         axis.title.x = element_text (size = 16, vjust = 0.2),
#         legend.position = "none") +
#   geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.20)/6))) + 
#   geom_hline(yintercept = -log10(Pthresh), linetype = "dashed") + 
#   geom_hline(yintercept = -log10(Pthresh), linetype = "dashed") + 
#   labs(x = "Base Pair Position (Mb)", y = "-log10(P)") +
#   scale_x_continuous(breaks = seq(0, 120e6, 2.5e6), labels = seq(0, 120, 2.5)) +
#   scale_y_continuous(breaks = seq(0, 12, 2), labels = seq(0, 12, 2))





snplistvec.150 <- read.table(paste0("results/2_Regional_Heritability_Run_a.ps_150SNP_ChrChrall_trans.txt"), header = T)
snplistvec.150$Chr <- factor(snplistvec.150$Chr)
snplistvec.150$Colour <- factor(snplistvec.150$Colour)
snplistvec.150$Label <- recoderFunc(snplistvec.150$Sex,
                                    c("Both", "Female", "Male"),
                                    c("a. Both Sexes", "b. Female", "c. Male"))


pdf("../Recombination Rate Manuscript/150snpwindowsRh2trans.pdf", width = 8, height = 8)
ggplot(snplistvec.150, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2, aes(shape = Error)) +
  theme_bw() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.150)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()

arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6)), Reg.h2.constraint == "Positive"), Sex, MedianCumuPos)
arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6)), Chr == 1), Sex)


snplistvec.50 <- read.table(paste0("results/2_Regional_Heritability_Run_b.ps_50SNP_ChrChrall_trans.txt"), header = T)
snplistvec.50$Chr <- factor(snplistvec.50$Chr)
snplistvec.50$Colour <- factor(snplistvec.50$Colour)
snplistvec.50$Label <- recoderFunc(snplistvec.50$Sex,
                                   c("Both", "Female", "Male"),
                                   c("a. Both Sexes", "b. Female", "c. Male"))


pdf("../Recombination Rate Manuscript/50snpwindowsRh2trans.pdf", width = 8, height = 8)
ggplot(snplistvec.50, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2, aes(shape = Error)) +
  theme_bw() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.50)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()


snplistvec.20 <- read.table(paste0("results/2_Regional_Heritability_Run_c.ps_20SNP_ChrChrall_trans.txt"), header = T)
snplistvec.20$Chr <- factor(snplistvec.20$Chr)
snplistvec.20$Colour <- factor(snplistvec.20$Colour)
snplistvec.20$Label <- recoderFunc(snplistvec.20$Sex,
                                   c("Both", "Female", "Male"),
                                   c("A. Both Sexes", "B. Females", "C. Males"))

pdf("../Recombination Rate Manuscript/20snpwindowsRh2trans.pdf", width = 8, height = 8, useDingbats = F)
ggplot(snplistvec.20, aes(MedianCumuPos, -log10(P), col = Colour, group = Chr)) + 
  geom_point(size = 2) +
  #geom_point(size = 2, aes(shape = Error)) +
  #geom_line(data = snplistvec.150, aes(MedianCumuPos, -log10(P)), alpha = 0.5, colour = "grey") + 
  theme_bw() +
  geom_line() +
  scale_color_manual(values = c("grey20", "grey55")) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
        axis.title.y = element_text (size = 14, vjust = 0.9),
        axis.title.x = element_text (size = 14, vjust = 0.2),
        strip.background = element_blank()) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05/(nrow(snplistvec.20)/6))) + 
  scale_x_continuous(breaks = tapply(maptab$Cumu, maptab$Chr, median)[1:26],
                     labels = c(1:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", 23, "", "", 26)) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)
dev.off()

head(snplistvec.20)

res150 <- arrange(filter(snplistvec.150, P < (0.05/(nrow(snplistvec.150)/6))), Sex)
res50 <- arrange(filter(snplistvec.50, P < (0.05/(nrow(snplistvec.50)/6))), Sex)
res20 <- arrange(filter(snplistvec.20, P < (0.05/(nrow(snplistvec.20)/6))), Sex)


filter(snplistvec.20,  V1 == "chr6.c.ps.1827.20")


regh2.results <- rbind(res150, res50, res20)
regh2.results <- filter(regh2.results, snpcount %in% c(20, 50, 150))
regh2.results <- select(regh2.results, V1, Sex, Chr, FirstPos, LastPos, Reg.h2, Reg.h2.SE, Genome.Rest.h2, Genome.Rest.h2.SE, P, Error)
filter(regh2.results, Reg.h2 > 0.01)


fullregh2 <- rbind(cbind(Window = 150, snplistvec.150),
                   cbind(Window = 50, snplistvec.50),
                   cbind(Window = 20, snplistvec.20))

head(fullregh2)
fullregh2 <- select(fullregh2,Window, Sex,  FirstPos, LastPos, Chr,  
                    Genome.Rest.h2, Genome.Rest.h2.SE, Reg.h2, Reg.h2.SE, 
                    Error, Reg.h2.constraint, Genome.GRM.LL, Genome.Reg.GRM.LL, Chi.Value, P)

names(fullregh2) <- c("Window.Size", "Sex", "First.SNP.Pos", "Last.SNP.Pos", "Chromosome",
                      "RestofGenome.GRM.h2", "RestofGenome.GRM.h2.SE", 
                      "Window.GRM.h2", "Window.GRM.h2.SE", "ASReml.Error",
                      "Window.GRM.h2.constraint", "RestofGenome.GRM.LL", 
                      "RestofGenome.and.Window.GRM.LL", "Chi.Value", "P")



write.table(fullregh2, "../Recombination Rate Manuscript/SupplementaryInformation/RegionalHeritabilityResultsTrans.csv", row.names = F, sep = ",", quote = F) 

