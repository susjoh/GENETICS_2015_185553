
library(reshape)
library(plyr)
library(ggplot2)

setwd("mixedmodelGWAS")
load("GLMM.data.Rdata")

source("C:/Users/Susan Johnston/Desktop/R Functions/recoderFunc.R")
source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")

waldfiles <- dir()[grep("Wald", dir())]
#waldfiles <- waldfiles[-grep("trans", waldfiles)]

effectfiles <- dir()[grep("EffectRes", dir())]
#effectfiles <- effectfiles[-grep("trans", effectfiles)]



waldresults <- NULL

for(i in waldfiles){
  waldtemp <- read.table(i, header = T)
  if(length(grep("trans", i)) > 0){
    waldtemp <- cbind(waldtemp, Type = "trans")
  } else {
    waldtemp <- cbind(waldtemp, Type = "cistrans")
  }
  waldresults <- rbind(waldresults, waldtemp)
  rm(waldtemp)
}

effectresults <- NULL

for(i in effectfiles){
  effecttemp <- read.table(i, header = T)
  if(length(grep("trans", i)) > 0){
    effecttemp <- cbind(effecttemp, Type = "trans")
  } else {
    effecttemp <- cbind(effecttemp, Type = "cistrans")
  }
  effectresults <- rbind(effectresults, effecttemp)
  rm(effecttemp)
}


write.table(waldresults, "waldresults.unprocessed.txt", row.names = F, sep = "\t", quote = F)
write.table(effectresults, "effectresults.unprocessed.txt", row.names = F, sep = "\t", quote = F)


maptab <- read.table("../results/2_Merged_map_g.txt", header = T)
maptab <- arrange(maptab, Chr, GenomicPosition)
maptab$Cumu <- cumsum(c(maptab$GenomicPosition[1], as.numeric(maptab$GenomicDiff)))[1:(nrow(maptab))]

maptab$Cumu2 <- maptab$Cumu + (25000000 * (maptab$Chr-1))


waldstore <- waldresults
waldresults <- waldstore

tail(waldresults)

wald.all.cistrans <- droplevels(subset(waldresults, Model == "All" & Type == "cistrans"))
wald.m.cistrans <- droplevels(subset(waldresults, Model == "Male" & Type == "cistrans"))
wald.f.cistrans <-  droplevels(subset(waldresults, Model == "Female" & Type == "cistrans"))

wald.all.trans <- droplevels(subset(waldresults, Model == "All" & Type == "trans"))
wald.m.trans <- droplevels(subset(waldresults, Model == "Male" & Type == "trans"))
wald.f.trans <-  droplevels(subset(waldresults, Model == "Female" & Type == "trans"))



lambda.all.cistrans <- median(wald.all.cistrans$Wald.statistic)/qchisq(0.5,2)
lambda.m.cistrans <- median(wald.m.cistrans$Wald.statistic)/qchisq(0.5,2)
lambda.f.cistrans <- median(wald.f.cistrans$Wald.statistic)/qchisq(0.5,2)

lambda.all.trans <- median(wald.all.trans$Wald.statistic)/qchisq(0.5,2)
lambda.m.trans <- median(wald.m.trans$Wald.statistic)/qchisq(0.5,2)
lambda.f.trans <- median(wald.f.trans$Wald.statistic)/qchisq(0.5,2)

waldresults$Label <- recoderFunc(waldresults$Model,
                                   c("All", "Female", "Male"),
                                   c("A. Both Sexes", "B. Females", "C. Males"))



waldresults <- join(waldresults, maptab)
waldresults <- subset(waldresults, Chr!= 0)

Pthresh <- 0.05/22273.61

chrinfo <- NULL

for(i in unique(sort(as.numeric(maptab$Chr)))){
  
  temp1 <- subset(maptab, Chr == i)
  temp2 <- data.frame(Chr = i,
                      Start = temp1[1,"Cumu2"],
                      Stop = temp1[nrow(temp1),"Cumu2"])
  chrinfo <- rbind(chrinfo, temp2)
}

chrinfo <- arrange(chrinfo, as.numeric(Chr))
chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$Chromosome2 <- c(1:10,"", 12, "", 14, "", 16, "", 18, "", "", 21, "",  "", 24, "", "", "X")

waldresults <- subset(waldresults, Chr!= 0)


waldresults$Colour <- "1"
waldresults$Colour[which(waldresults$Chr %in% seq(0, 26, 2))] <- "2"

# Plot

gwas <- ggplot(subset(waldresults, Type = "cistrans"), aes(Cumu2, -log10(Pr.Chisq.), col = factor(Colour))) +
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
                     labels = chrinfo$Chromosome2) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Chromosome", y = "-log10(P)") +
  facet_wrap(~Label, ncol = 1)

gwa_res.null<-NULL

for(i in sort(unique(waldresults$Label))){
  subdata <- subset(waldresults, Label == i)
  x <- data.frame(obs = sort(-log10(subdata$Pr.Chisq.)),
                  exp = sort(-log10(pchisq(qchisq(seq(1/nrow(subdata),1,1/nrow(subdata)), 2), 2))),
                  Label = i)
  gwa_res.null <- rbind(gwa_res.null, x)
  rm(subdata, x)
}
gc()


ppplot <- 
  ggplot(gwa_res.null, aes(exp, obs)) +
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
  #coord_cartesian(xlim = c(0, 30)) +
  scale_y_continuous(breaks = seq(0, 30, 2)) +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
  facet_wrap(~Label, ncol = 1)


pdf("../../Recombination Rate Manuscript/GWAScistransv2.pdf", width = 10, height = 8, useDingbats = F)
multiplot(gwas, ppplot, cols = 2, layout = matrix(c(1, 1, 2), nrow = 1))
dev.off()


#~~ Make results table

waldresults <- arrange(waldresults, Pr.Chisq.)
waldresults[which(waldresults$Pr.Chisq. < Pthresh),]

waldresults <- dplyr::select(waldresults, Model, Type, SNP.Name, Chr, GenomicPosition, Df, Wald.statistic, Pr.Chisq.)
head(waldresults)


gwas.cistrans <- read.table("../results/2_GWAS_Results_cistrans.txt", header = T)
head(gwas.cistrans)
gwas.cistrans$Model <- recoderFunc(gwas.cistrans$Analysis,
                                   c("BothSexes", "Females", "Males"),
                                   c("All", "Female", "Male"))

gwas.cistrans <- subset(gwas.cistrans, select = c(Model, A1, A2, N, SNP.Name, CallRate, MAF))
waldresults <- join(waldresults, gwas.cistrans)
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

effectresults2 <- melt(effectresults, id.vars = c("SNP.Name", "Model", "Effect2", "Type"))
head(effectresults2)
effectresults2$NewEffect <- paste(effectresults2$Effect2, effectresults2$variable, sep = ".")
effectresults2 <- subset(effectresults2, select = -c(Effect2, variable))

effectresults3 <- cast(effectresults2, SNP.Name + Model + Type ~ NewEffect)

head(effectresults3)

head(waldresults)

waldresults <- join(waldresults, effectresults3)
waldresults <- arrange(waldresults, Type, Model, Chr, GenomicPosition)

write.table(waldresults, "../../Recombination Rate Manuscript/SupplementaryInformation/GWASresultsASREML.csv", row.names = F, quote = F, sep = ",")

