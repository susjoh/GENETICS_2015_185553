```{r set-options, echo=FALSE, cache=FALSE}
options(width=120)
options(error = function(){    # Beep on error
  beepr::beep()
  Sys.sleep(1)
  }
 )

.Last <- function() {          # Beep on exiting session
  beepr::beep()
  Sys.sleep(1)
  }
```
 

```{r}
setwd("simulationsRR/")

library(reshape2)
library(plyr)
library(ggplot2)

corr.res <- data.frame(Iteration = c(1:100),
                       AutosomeSlope = NA,
                       AutosomeAdjR2 = NA,
                       TotalAutosomeSlope = NA,
                       TotalAutosomeAdjR2 = NA)

rectab <- read.table("../results/2_IndivRR_FullClean_g.txt", header = T)


for(i in corr.res$Iteration){
  load(paste0("recsummfiles/SummaryModels", i, ".Rdata"))
  rownumber <- which(corr.res$Iteration == i)
  corr.res$AutosomeSlope[rownumber] <- summary(autosomeCorr)$coefficients[2,1]
  corr.res$AutosomeAdjR2[rownumber] <- summary(autosomeCorr)$adj.r.squared
  
  fit <- lm(recsumm$TrueRecombCount ~ recsumm$TotalRecombCount)
  
  corr.res$TotalAutosomeSlope[rownumber] <- summary(fit)$coefficients[2, 1]
  corr.res$TotalAutosomeAdjR2[rownumber] <- summary(fit)$adj.r.squared
  
  
  head(recsumm)
  
  rm(allchrCorr, autosomeCorr, xchrCorr, recsumm, fit)
  gc()
}

# corr.res
corr.res2 <- melt(corr.res, id.vars = "Iteration")

head(corr.res2)


corr.res2$Label[which(corr.res2$variable == "AutosomeSlope")] <- "Chromtids: Slope"
corr.res2$Label[which(corr.res2$variable == "AutosomeAdjR2")] <-  "Chromatids: Adj. R2"
corr.res2$Label[which(corr.res2$variable == "TotalAutosomeSlope")] <- "Individuals: Slope"
corr.res2$Label[which(corr.res2$variable == "TotalAutosomeAdjR2")] <-  "Individuals: Adj. R2"

corr.summ <- data.frame(Mean = tapply(corr.res2$value, corr.res2$variable, mean),
                        SD   = tapply(corr.res2$value, corr.res2$variable, sd))
corr.summ$variable <- row.names(corr.summ)

corr.res2 <- join(corr.res2, corr.summ)

corr.summ

require(grid)
ggplot(corr.res2, aes(value)) + 
  geom_histogram(col = "black", fill = "grey") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 10, vjust = 0),
        axis.text.y  = element_text (size = 10),
        strip.text.x = element_text (size = 12, vjust = 0.7),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12, vjust = 0.2),
        panel.margin = unit(2, "lines"),
        strip.background = element_blank()) +
  facet_wrap(~Label, scales = "free_x") +
  geom_vline(aes(xintercept = Mean), linetype = "dashed")


#~~~~~~~~~~~~~~~~
# indiv analysis
#~~~~~~~~~~~~~~~~

recsumm.list <- list()

for(i in corr.res$Iteration){
  # print(i)
  load(paste0("recsummfiles/SummaryModels", i, ".Rdata"))
  recsumm.list[[i]] <- recsumm
  rm(allchrCorr, autosomeCorr, xchrCorr, recsumm)
  gc()
}



id.frame <- data.frame(RRID.Offspring = recsumm.list[[1]]$RRID,
                       AdjR2 = NA,
                       Slope = NA,
                       MeanInfLoci = NA,
                       SDInfLoci = NA)
id.frame$RRID.Offspring <- as.character(id.frame$RRID.Offspring)

head(id.frame)


for(i in 1:nrow(id.frame)){
  if(i %in% seq(1, nrow(id.frame), 100)) print(paste("Analysing row", i, "of", nrow(id.frame)))
  
  tempframe <- data.frame(Iteration = corr.res$Iteration,
                          TotalRecombCount = NA,
                          TrueRecombCount = NA,
                          TotalInfLoci = NA)
  
  for(j in corr.res$Iteration){
    rownumber <- which(tempframe$Iteration == j)
    tempframe$TotalRecombCount[rownumber] <- recsumm.list[[j]][which(recsumm.list[[j]]$RRID == id.frame$RRID.Offspring[i]),"TotalRecombCount"][[1]]
    tempframe$TrueRecombCount[rownumber]  <- recsumm.list[[j]][which(recsumm.list[[j]]$RRID == id.frame$RRID.Offspring[i]),"TrueRecombCount" ][[1]] 
    tempframe$TotalInfLoci[rownumber]  <- recsumm.list[[j]][which(recsumm.list[[j]]$RRID == id.frame$RRID.Offspring[i]),"TotalInfLoci" ][[1]] 
    
    
    rm(rownumber)
  }
  
  fit <- lm(TotalRecombCount ~ TrueRecombCount, tempframe)
  id.frame$AdjR2[i] <- summary(fit)$adj.r.squared
  id.frame$Slope[i] <- summary(fit)$coefficients[2,1]
  id.frame$MeanInfLoci[i] <- mean(tempframe$TotalInfLoci)
  id.frame$SDInfLoci[i] <- sd(tempframe$TotalInfLoci)
  
}

head(id.frame)

ggplot(id.frame, aes(AdjR2)) + geom_histogram()
ggplot(id.frame, aes(Slope)) + geom_histogram()

head(id.frame)
id.frame$RRID <- unlist(sapply(id.frame$RRID.Offspring, function(x) strsplit(x, split = "_")[[1]][1]))
id.frame$Offspring.ID <- unlist(sapply(id.frame$RRID.Offspring, function(x) strsplit(x, split = "_")[[1]][2]))

head(recsumm.list[[1]])
test <- read.table("../results/2_TotalSummIndivRR_FullClean_g.txt", header = T)
head(test)
pedigree <- read.table("../data/pedigree_20130920.txt", header = T)
names(pedigree)[1] <- "RRID"

#id.frame <- join(id.frame, 

inbreeding <- read.table("../gcta/20150129_autosomal_IBD.ibc", header = T)
inbreeding <- inbreeding[,c("IID", "Fhat3")]
names(inbreeding)[1] <- "RRID"

id.frame <- join(id.frame, test)
head(id.frame)

ggplot(id.frame, aes(MeanInfLoci)) + geom_histogram()
ggplot(id.frame, aes(MeanInfLoci, AdjR2)) + geom_point()
ggplot(id.frame, aes(MeanInfLoci, TotalInfLoci)) + geom_point()
ggplot(id.frame, aes(MeanInfLoci, TotalInfLoci)) + geom_point()
ggplot(id.frame, aes(RRID.Fhat3, MeanInfLoci)) + geom_point()
ggplot(id.frame, aes(Offspring.Fhat3, AdjR2)) + geom_point()
ggplot(id.frame, aes(RRID.Fhat3, AdjR2)) + geom_point()

lifefit <- read.table("../results/2_Lifetime_Fitness_Measures_20150313.txt", header = T, sep = "\t")
head(lifefit)
lifefit <- subset(lifefit, select = c(ID, TotalOffspringBorn))
names(lifefit)[1] <- "RRID"

id.frame <- join(id.frame, lifefit)
ggplot(id.frame, aes(TotalOffspringBorn, AdjR2)) + geom_point()

# id.frame[which(id.frame$AdjR2 < 0.95),]
rm(recsumm.list)
gc()

# rectab[which(rectab$RRID == 7139),"data.v2"]

write.table(id.frame[,c("RRID", "UniqueID2", "AdjR2", "Slope", "MeanInfLoci", "SDInfLoci")],
            "../results/4_SimulationResults.txt", row.names = F, sep = "\t", quote = F)
names(id.frame)

mean(id.frame$AdjR2)
```

```{r}

#~~ set analysis parameters

AnalysisSuffix <- "g"

#~~ read in data

# Individual summed data
recsumm  <- read.table(paste0("../results/2_TotalSummIndivRR_FullClean_", AnalysisSuffix, ".txt"), header = T, sep = "\t")

# Pedigree
pedigree <- read.table("../data/pedigree_20130920.txt", header = T)

# Simulation Results

simtab <- read.table("../results/4_SimulationResults.txt", header = T)

basedata <- read.table("../data/20150303_SoayBaseData.txt", header = T, sep = "\t")
inbreeding <- read.table("../gcta/20150129_autosomal_IBD.ibc", header = T)

head(basedata)
names(basedata)[1] <- "RRID"
basedata <- subset(basedata, Sex %in% c(1, 2))
basedata$SEX <- ifelse(basedata$Sex == 1, "Female", "Male")
head(simtab)
head(inbreeding)
names(inbreeding)[2] <- "RRID"

require(plyr)
require(ggplot2)

simtab <- join(simtab, basedata[,c("RRID", "SEX")])
simtab <- join(simtab, inbreeding[,c("RRID", "Fhat3")])
source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")
multiplot(
ggplot(simtab, aes(Fhat3, AdjR2, col = SEX)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal") +
  labs(x = expression(paste("Genomic Inbreeding Coefficient (", hat("F"), ")")), y = expression("Adjusted R"^" 2")),


ggplot(simtab, aes(Fhat3, Slope, col = SEX)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal") +
  labs(x = expression(paste("Genomic Inbreeding Coefficient (", hat("F"), ")")), y = "Linear Regression Slope"),
cols = 2)

ggplot(simtab, aes(SEX, AdjR2, col = SEX)) + geom_boxplot(notch = T) + theme_bw() +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 strip.background = element_blank(),
                legend.position = "top",
                legend.direction = "horizontal") +
    scale_colour_brewer(palette = "Set1") +
  labs(x = "Sex", y = expression("Adjusted R"^" 2"))

fit <- glm(AdjR2 ~ SEX + Fhat3, data = simtab)
summary(fit)
fit <- glm(Slope ~ SEX + Fhat3, data = simtab)
summary(fit)
require(MASS)
library(hglm)
fit <- hglm2(fixed = simtab$AdjR2 ~ simtab$Fhat3 + simtab$SEX, random = ~ 1|simtab$RRID, family = Beta(link = "logit"))
summary(fit)

```