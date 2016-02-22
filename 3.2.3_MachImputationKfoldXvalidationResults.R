# Analyse output from cross-validation
# Author: Susan Johnston


library(beepr)
library(ggplot2)

simtab <- NULL

xval <- 10

tenfoldvec <- dir("mach/xvalidation")[grep("10fold", dir("mach/xvalidation"))]
tenfoldvec <- gsub("simulation_10fold_Run", "", tenfoldvec)
tenfoldvec <- gsub(".Rdata", "", tenfoldvec)
tenfoldvec <- sort(as.numeric(tenfoldvec))

valid.list <- list()

# for(run in tenfoldvec){
# 
#   if(run %in% seq(1, 1000, 10)) print(paste("Analysing run", run)) 
#   
#   load(paste0("mach/xvalidation/simulation_", xval, "fold_Run", run, ".Rdata"))
#   valid.list[[run]] <- valid.ids
#   simtab <- rbind(simtab, cbind(Iteration = run, restab.subset))
#   
#   rm(i, id.vec, imputed, imputed.haplos, imputed.pr, imputed.snplist, 
#      j, mapfile, pedfile, pedfile.subset, restab.subset,
#      tabtomelt, valid.ids, valid.tab, xvalidation)
# }

load("test.Rdata")

head(simtab)
#simtab <- subset(simtab, select = -Link)
gc()
beep()
ggplot(simtab, aes(Dosage)) + geom_histogram(binwidth = 0.05)

simtab2 <- simtab
simtab <- simtab2

length(unique(simtab$id))

#~~ Rescore Dosage between 0.25 and 0.75, 1.25 and 1.75

simtab$Genotype[which(simtab$Dosage > 0.25 & simtab$Dosage < 0.75)] <- "0/0"
simtab$Genotype[which(simtab$Dosage > 1.25 & simtab$Dosage < 1.75)] <- "0/0"

head(simtab)

table(simtab$Genotype, simtab$value)

#~~ Per locus

simtab$Match <- ifelse(simtab$Genotype == simtab$value, 1, 0)
simtab2 <- subset(simtab, Genotype != "0/0")

locus.res <- data.frame(Matches = tapply(simtab$Match, simtab$SNP.Name, sum),
                        Count = tapply(simtab$Match, simtab$SNP.Name, length),
                        CountNonMissing = tapply(simtab2$Match, simtab2$SNP.Name, length))
locus.res$ProportionMatch <- locus.res$Matches/locus.res$CountNonMissing
locus.res$ProportionScored <- locus.res$CountNonMissing/locus.res$Count
locus.res$SNP.Name <- row.names(locus.res)

head(locus.res)
ggplot(locus.res, aes(ProportionMatch)) + geom_histogram()
ggplot(locus.res, aes(ProportionScored)) + geom_histogram()
ggplot(locus.res, aes(ProportionMatch, ProportionScored)) + geom_point()


write.table(locus.res, "mach/10foldXvalidationAccuracy.Locus.txt", sep = "\t", quote = F, row.names = F)
# 
# simtab <- subset(simtab, SNP.Name %in% locus.res$SNP.Name[which(locus.res$ProportionMatch > 0.99)])
# simtab2 <- subset(simtab, Genotype != "0/0")


#~~ Per ID

id.res <- data.frame(Matches = tapply(simtab$Match, simtab$id, sum),
                      Count = tapply(simtab$Match, simtab$id, length),
                      CountNonMissing = tapply(simtab2$Match, simtab2$id, length))
id.res$ProportionMatch <- id.res$Matches/id.res$CountNonMissing
id.res$ProportionScored<- id.res$Matches/id.res$Count
id.res$id <- row.names(id.res)

head(id.res)
ggplot(id.res, aes(Count)) + geom_histogram()
ggplot(id.res, aes(CountNonMissing)) + geom_histogram()
ggplot(id.res, aes(ProportionMatch)) + geom_histogram()
ggplot(id.res, aes(ProportionScored)) + geom_histogram()

plot(id.res$ProportionScored, id.res$ProportionMatch)

write.table(id.res, "mach/10foldXvalidationAccuracy.ID.txt", sep = "\t", quote = F, row.names = F)


#~~ Per ID per simulation
library(reshape)
id.sim.res1 <- data.frame(Matches = tapply(simtab$Match,
                                           list(simtab$id, simtab$Iteration),
                                           sum))
id.sim.res1$id <- row.names(id.sim.res1)
id.sim.res1 <- melt(id.sim.res1, id.vars = "id")
id.sim.res1 <- subset(id.sim.res1, !is.na(value))
head(id.sim.res1)
names(id.sim.res1) <- c("id", "Iteration", "Matches")
id.sim.res1$Iteration <- gsub("Matches.", "", id.sim.res1$Iteration)

id.sim.res2 <- data.frame(Count = tapply(simtab$Match,
                                         list(simtab$id, simtab$Iteration),
                                         length))
id.sim.res2$id <- row.names(id.sim.res2)
id.sim.res2 <- melt(id.sim.res2, id.vars = "id")
id.sim.res2 <- subset(id.sim.res2, !is.na(value))
head(id.sim.res2)
names(id.sim.res2) <- c("id", "Iteration", "Count")
id.sim.res2$Iteration <- gsub("Count.", "", id.sim.res2$Iteration)


id.sim.res3 <- data.frame(Count = tapply(simtab2$Match,
                                         list(simtab2$id, simtab2$Iteration),
                                         length))
id.sim.res3$id <- row.names(id.sim.res3)
id.sim.res3 <- melt(id.sim.res3, id.vars = "id")
id.sim.res3 <- subset(id.sim.res3, !is.na(value))
head(id.sim.res3)
names(id.sim.res3) <- c("id", "Iteration", "CountNonMissing")
id.sim.res3$Iteration <- gsub("Count.", "", id.sim.res3$Iteration)

library(plyr)
id.sim.res <- join(id.sim.res1, id.sim.res2)
id.sim.res <- join(id.sim.res , id.sim.res3)
head(id.sim.res)

rm(id.sim.res1, id.sim.res2, id.sim.res3)

id.sim.res[1:100,]
id.sim.res$ProportionMatch <- id.sim.res$Matches/id.sim.res$CountNonMissing
id.sim.res$ProportionScored <- id.sim.res$CountNonMissing/id.sim.res$Count

head(id.sim.res)

gc()

ggplot(id.sim.res, aes(ProportionScored)) + geom_histogram()
ggplot(id.sim.res, aes(ProportionMatch)) + geom_histogram()

ggplot(id.sim.res, aes(ProportionScored, ProportionMatch)) + geom_point(alpha = 0.05)

ggplot(id.sim.res, aes(id, ProportionMatch)) + geom_point(alpha = 0.2)

id.sim.res.summ <- data.frame(MeanPropScored = tapply(id.sim.res$ProportionScored, id.sim.res$id, mean),
                              MeanPropMatch  = tapply(id.sim.res$ProportionMatch, id.sim.res$id, mean))

ggplot(id.sim.res.summ, aes(MeanPropScored, MeanPropMatch)) + geom_point(alpha = 0.2)
# save.image("test.Rdata")

id.sim.res.summ[which(id.sim.res.summ$MeanPropMatch < 0.95),]

# remove ID's where proportion of loci scored was less than 0.95

pedigree <- read.table("data/pedigree_20130920.txt", header = T)

pedigree[which(pedigree$ANIMAL == row.names(id.sim.res.summ[which(id.sim.res.summ$MeanPropMatch < 0.95),]
)),]

pedigree[which(pedigree$ANIMAL %in% row.names(id.sim.res.summ)),]



simtab <- subset(simtab, id %in% id.res$id[which(id.res$ProportionScored > 0.95)])

#~~ Redo the per locus information

simtab$Match <- ifelse(simtab$Genotype == simtab$value, 1, 0)
simtab2 <- subset(simtab, Genotype != "0/0")

locus.res <- data.frame(Matches = tapply(simtab$Match, simtab$SNP.Name, sum),
                        Count = tapply(simtab$Match, simtab$SNP.Name, length),
                        CountNonMissing = tapply(simtab2$Match, simtab2$SNP.Name, length))
locus.res$ProportionMatch <- locus.res$Matches/locus.res$CountNonMissing
locus.res$ProportionScored <- locus.res$CountNonMissing/locus.res$Count
locus.res$SNP.Name <- row.names(locus.res)

head(locus.res)
ggplot(locus.res, aes(ProportionMatch)) + geom_histogram()
ggplot(locus.res, aes(ProportionScored)) + geom_histogram()











