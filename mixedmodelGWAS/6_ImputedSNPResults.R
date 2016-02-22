
library(reshape)
library(plyr)
library(ggplot2)
library(GenABEL)

setwd("mixedmodelGWAS")

source("C:/Users/Susan Johnston/Desktop/R Functions/recoderFunc.R")
source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")

imputedIDs <- read.table("../mach/chr6HDsegPostQC.ped")$V2


impabel <- load.gwaa.data(genofile = "../mach/imputedPLINK.genabel",
                          phenofile = "../mach/imputedPLINK.phenofile")
impabel.original <- impabel[which(idnames(impabel) %in% imputedIDs),]

waldresults <- read.table("ImputedWaldtest.It_1_126_cistrans.txt", header = T)
effectresults <- read.table("ImputedEffectRes.It_1_126_cistrans.txt", header = T)


waldresults$Label <- recoderFunc(waldresults$Model,
                                 c("All", "Female", "Male"),
                                 c("A. Both Sexes", "B. Females", "C. Males"))

chr6results <- read.table("../mach/ImputationDataForPlotting.txt", header = T)
head(chr6results)

chr6results <- subset(chr6results, select = c(Position, Array, SNP.Name))
waldresults <- join(chr6results, waldresults)
head(waldresults)

Pthresh <- 0.05/22273.61
waldresults$Pr.Chisq.2 <- pchisq(waldresults$Wald.statistic, waldresults$Df, lower.tail = F)

effectresults[which(effectresults$SNP.Name == "oar3_OAR6_116402578"),]


ggplot() + 
  geom_point(data = subset(waldresults, Model == "Female"), aes(Position, -log10(Pr.Chisq.2), colour = Array, shape = Array), size = 2) +
  scale_colour_manual(values = c("black", "red")) +
  geom_point(data = subset(waldresults, Array == "snp50" & Model == "Female"), aes(Position, -log10(Pr.Chisq.2)), size = 2, col = "black") +
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





#~~ Make results table


waldresults <- arrange(waldresults, Pr.Chisq.)
waldresults[which(waldresults$Pr.Chisq. < Pthresh),]

waldresults <- dplyr::select(waldresults, Model, Type, SNP.Name, Chr, GenomicPosition, Df, Wald.statistic, Pr.Chisq.)
head(waldresults)

snp.info <- summary.snp.data(gtdata(impabel))
head(snp.info)
snp.info$SNP.Name <- row.names(snp.info)

snp.info <- subset(snp.info, select = c(A1, A2, SNP.Name, CallRate, Q.2))

snp.info2  <- summary.snp.data(gtdata(impabel.original))
snp.info2$SNP.Name <- row.names(snp.info2)

head(snp.info2)
snp.info2 <- subset(snp.info2, select = c(SNP.Name, Q.2))
names(snp.info2) <- c("SNP.Name", "MAF.Original")

snp.info <- join(snp.info, snp.info2)

waldresults <- join(waldresults, snp.info)
head(waldresults)



head(effectresults)
effectresults <- subset(effectresults, Effect != "(Intercept)")
effectresults$Effect <- gsub("as.factor\\(SNP\\)_", "", effectresults$Effect)
effectresults <- join(effectresults, unique(waldresults[,c("SNP.Name", "A1", "A2")]))
effectresults <- subset(effectresults, !is.na(A1))

# STILL TO DO

effectresults$temp1 <- paste0(effectresults$A1, "/", effectresults$A1)
effectresults$temp2 <- paste0(effectresults$A2, "/", effectresults$A2)

effectresults$Effect2 <- ifelse(effectresults$Effect == effectresults$temp1, "AA",
                                ifelse(effectresults$Effect == effectresults$temp2, "BB", "AB"))
head(effectresults)
effectresults <- subset(effectresults, select = -c(Effect, temp1, temp2, A1, A2))

effectresults2 <- melt(effectresults, id.vars = c("SNP.Name", "Model", "Effect2"))
head(effectresults2)
effectresults2$NewEffect <- paste(effectresults2$Effect2, effectresults2$variable, sep = ".")
effectresults2 <- subset(effectresults2, select = -c(Effect2, variable))

effectresults3 <- cast(effectresults2, SNP.Name + Model  ~ NewEffect)

head(effectresults3)

head(waldresults)

waldresults <- join(waldresults, effectresults3)
waldresults <- arrange(waldresults, Model, Position)

write.table(waldresults, "../../Recombination Rate Manuscript/SupplementaryInformation/GWASresultsASREMLImputed.csv", row.names = F, quote = F, sep = ",")



recsumm <- read.table("../results/2_TotalSummIndivRR_FullCleanPostSim_g.txt", header = T)

sigsnptab <- data.frame(SigGeno = as.character.gwaa.data(impabel[,c("oar3_OAR6_116402578")]))
sigsnptab$RRID <- row.names(sigsnptab)
recsumm <- join(recsumm, sigsnptab)
head(recsumm)

table(recsumm$RRID.SEX, recsumm$oar3_OAR6_116402578)

# table(recsumm$SigGeno.s74824.1, recsumm$SigGeno.s14138.1, useNA = "always")
# table(recsumm$SigGeno.s74824.1, recsumm$SigGeno.s10844.1, useNA = "always")
# table(recsumm$SigGeno.s14138.1, recsumm$SigGeno.s10844.1, useNA = "always")
# 


ggplot(recsumm, aes(oar3_OAR6_116402578, TotalRecombCount)) + 
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

geteffect <- read.table("ImputedEffectRes.It_1_126_cistrans.txt", header = T)
geteffect <- subset(geteffect, SNP.Name == "oar3_OAR6_116402578" & Model != "All")

geteffect2 <- geteffect[grep("Intercept", geteffect$Effect),]

geteffect <- geteffect[grep("SNP", geteffect$Effect),]
geteffect$Effect <- rep(c("A/A", "G/A", "G/G"), times = 2)
geteffect

geteffect$solution2 <- geteffect$solution
geteffect$solution2[which(geteffect$Effect == "A/A")] <- geteffect2$solution
geteffect$solution2se <- geteffect$std.error
geteffect$solution2se[which(geteffect$Effect == "A/A")] <- geteffect2$std.error

geteffect$solution3 <- NA
geteffect$solution3[1:3] <- geteffect$solution[1:3] + geteffect2$solution[1]
geteffect$solution3[4:6] <- geteffect$solution[4:6] + geteffect2$solution[2]

geteffect$Analysis <- rep(c("Males", "Females"), each = 3)

geteffect$Label <- rep(c("B. Males", "A. Females"), each = 3)

temp <- data.frame(table(recsumm$RRID.SEX, recsumm$oar3_OAR6_116402578))
names(temp)[3] <- "Observations"

temp2 <- unique(subset(recsumm, select = c(RRID, RRID.SEX, oar3_OAR6_116402578)))

temp3 <- data.frame(table(temp2$RRID.SEX, temp2$oar3_OAR6_116402578))
names(temp3)[3] <- "UniqueIDs"

temp <- join(temp, temp3)
names(temp)[1:2] <- c("Analysis", "Effect")
temp$Analysis <- paste0(temp$Analysis, "s")
temp$Effect <- as.character(temp$Effect)

geteffect <- join(geteffect, temp)
rm(temp, temp2, temp3)


geteffect.m <- subset(geteffect, Analysis == "Males")
geteffect.f <- subset(geteffect, Analysis == "Females")

pdf("../../Recombination Rate Manuscript/snpeffectsizes.pdf", width = 6, height = 4, useDingbats = F)
#postscript("../Recombination Rate Manuscript/snpeffectsizes.ps", width = 600, height = 400)

multiplot(
  ggplot(geteffect.f, aes(Effect, solution3)) + 
    geom_rect(xmin = 0, xmax = 4, ymin = max(geteffect.f$solution3) + 1.8, ymax = max(geteffect.f$solution3) + 3, fill = "white", col = "black") +
    geom_point(size = 4) + 
    geom_errorbar(aes(ymax = solution3 + solution2se, ymin = solution3 - solution2se), width = 0) +
    facet_wrap(~Label) +
    labs(x = "oar3_OAR6_116402578", y = "Effect size") +
    theme_bw() +
    coord_cartesian(ylim = c(max(geteffect.f$solution3) - 7, max(geteffect.f$solution3) +3)) +
    geom_text(aes(y = max(geteffect.f$solution3) + 2.7, label = Observations), size = 4) +
    geom_text(aes(y = max(geteffect.f$solution3) + 2.2, label = paste0("(", UniqueIDs, ")")), size = 4) +
    scale_y_continuous(breaks = seq(0, 29, 1)) +
    theme(axis.text.x  = element_text (size = 12),
          axis.text.y  = element_text (size = 12),
          strip.text.x = element_text (size = 14, vjust = 0.7, hjust = 0),
          axis.title.y = element_text (size = 14, vjust = 1.4),
          axis.title.x = element_text (size = 14, vjust = 0.2),
          strip.background = element_blank())
  ,
  ggplot(geteffect.m, aes(Effect, solution3)) + 
    geom_rect(xmin = 0, xmax = 4, ymin = max(geteffect.m$solution3) + 1.8, ymax = max(geteffect.m$solution3) + 3, fill = "white", col = "black") +
    geom_point(size = 4) + 
    geom_errorbar(aes(ymax = solution3 + solution2se, ymin = solution3 - solution2se), width = 0) +
    facet_wrap(~Label) +
    labs(x = "oar3_OAR6_116402578", y = "Effect size") +
    theme_bw() +
    coord_cartesian(ylim = c(max(geteffect.m$solution3) - 7, max(geteffect.m$solution3) +3)) +
    geom_text(aes(y = max(geteffect.m$solution3) + 2.7, label = Observations), size = 4) +
    geom_text(aes(y = max(geteffect.m$solution3)+ 2.2, label = paste0("(", UniqueIDs, ")")), size = 4) +
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











