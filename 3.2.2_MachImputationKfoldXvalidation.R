# k-fold cross validation of imputed genotypes
# Author: Susan Johnston

i = 1
xvalidation <- 2

setwd("mach")

library(plyr)
library(dplyr)
library(reshape)

#~~ make dat file

mapfile <- read.table("testHD.dat", stringsAsFactors = F)

#~~ Sample IDs

pedfile <- read.table("testHD.ped")
id.vec <- pedfile$V2

valid.ids <- sample(id.vec, size = length(id.vec)/xvalidation, replace = F)



#~~ Create temporary file

write.table(pedfile[which(!pedfile$V2 %in% valid.ids),], paste0("temp", i, ".ped"), row.names = F, col.names = F, quote = F)

#~~ Run mach

system("cmd", input = paste0("mach1 -d testHD.dat -p temp", i, ".ped --rounds 20 --states 200 --phase --interim 5 --sample 5 --compact --prefix temp", i, "snpset"))
system("cmd", input = paste0("mach1 -d test50.dat -p test50.ped -s testImpute.snps -h temp", i, "snpset --crossover temp", i, "snpset.erate --errormap temp", i, "snpset.rec --greedy --rounds 10 --mle --mldetails --prefix temp", i, "imputed"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Read in imputed genotypes     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

imputed <- read.table(paste0("temp", i, "imputed.mlgeno"))[,-2]
imputed.pr <- read.table(paste0("temp", i, "imputed.mldose"))[,-2]
imputed.snplist <- read.table(paste0("temp", i, "imputed.mlinfo"), header = T, stringsAsFactors = F)
imputed.haplos <- read.table(paste0("temp", i, "snpset"))



names(imputed) <- c("Link", as.character(imputed.snplist$SNP))
names(imputed.pr) <- c("Link", as.character(imputed.snplist$SNP))

pr.melt <- melt(imputed.pr, id.vars = c("Link"))
head(pr.melt)
names(pr.melt)[2:3] <- c("SNP.Name", "Dosage")


geno.melt <- melt(imputed, id.vars = c("Link"))
head(geno.melt)
names(geno.melt)[2:3] <- c("SNP.Name", "Genotype")


restab <- join(geno.melt, pr.melt)

restab$id <- gsub("1->", "", restab$Link)

#~~ Get the old genotypes

restab.subset <- filter(restab, id %in% valid.ids)

pedfile.subset <- filter(pedfile, V2 %in% valid.ids)

tabtomelt <- data.frame(id = pedfile.subset[,2])

for(j in seq(6, ncol(pedfile.subset), 2)){
  tabtomelt <- cbind(tabtomelt, paste(pedfile.subset[,j], pedfile.subset[,j+1], sep = "/"))
}
names(tabtomelt)[2:ncol(tabtomelt)] <- mapfile$V2


valid.tab <- melt(tabtomelt, id.vars = "id")
names(valid.tab)[2] <- "SNP.Name"

restab.subset <- join(restab.subset, valid.tab)
table(restab.subset$Genotype, restab.subset$value)

restab.subset$Genotype <- as.character(restab.subset$Genotype)
restab.subset$value    <- as.character(restab.subset$value)


restab.subset$Genotype[which(restab.subset$Genotype ==  "C/A")] <- "A/C"
restab.subset$Genotype[which(restab.subset$Genotype ==  "G/A")] <- "A/G"
restab.subset$Genotype[which(restab.subset$Genotype ==  "T/A")] <- "A/T"
restab.subset$Genotype[which(restab.subset$Genotype ==  "G/C")] <- "C/G"
restab.subset$Genotype[which(restab.subset$Genotype ==  "T/C")] <- "C/T"
restab.subset$Genotype[which(restab.subset$Genotype ==  "T/G")] <- "G/T"

restab.subset$value[which(restab.subset$value ==  "C/A")] <- "A/C"
restab.subset$value[which(restab.subset$value ==  "G/A")] <- "A/G"
restab.subset$value[which(restab.subset$value ==  "T/A")] <- "A/T"
restab.subset$value[which(restab.subset$value ==  "G/C")] <- "C/G"
restab.subset$value[which(restab.subset$value ==  "T/C")] <- "C/T"
restab.subset$value[which(restab.subset$value ==  "T/G")] <- "G/T"

table(restab.subset$Genotype, restab.subset$value)

rm(geno.melt, pr.melt, restab)

save.image(paste0("simulation_", xvalidation, "fold_Run", i, ".Rdata"))

}