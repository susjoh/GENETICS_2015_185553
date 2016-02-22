#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   File to create working genetic files based on data updates  #
#   Susan Johnston January 2013                                 #
#   Updated after Database v1.80 beta March 2013                #
#   Updated with Oar3.1 positions April 2013                    #
#   Updated for analysis "g" January 2015                       #
#                                                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Requires plink, gawk, mkdir, rm in PATH

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment and load in data         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ Load libraries

library(ggplot2)
library(plyr)
library(GenABEL)


#~~ Load functions

source("0_Chrompic_Functions.R")
source("0_Crimap_Functions_Faster_v2.R")

#~~ read in data and format

basedata <- read.table("data/20140210_SoayBaseData.txt", header = T, sep = "\t")
pedigree <- read.table("data/pedigree_20130920.txt", header = T)
plink.file <- "data/20130410merged1_66nodups.oar3.1"
badsnps <- read.table("data/1_Wrongly_Mapped_SNPs.txt", header = T)


names(basedata) <- toupper(names(basedata))
names(pedigree) <- toupper(names(pedigree))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Conduct a sex check on the plink file               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("cmd", input = paste("plink --file", plink.file, " --check-sex --sheep --from s17709.1 --to OARX_95154058.1"))
sexcheck <- read.table("plink.sexcheck", header = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create GenABEL files                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

newsheep <- unique(basedata[,c("ID", "HORN", "SEX")])
names(newsheep) <- c("id", "NewHorn", "pedsex")
newsheep$pedsex <- recoderFunc(newsheep$pedsex, c(2, 1, 3), c(1, 0, NA))
newsheep$pedsex2 <- ifelse(!is.na(newsheep$pedsex), "known", "unknown")

#~~ now match to the sexcheck information

sexcheck$IID[which(sexcheck$IID[which(sexcheck$STATUS == "PROBLEM")] %in% pedigree$MOTHER)]

plot(sexcheck$F, col = sexcheck$PEDSEX)

sexcheck$Is.Mum <- ifelse(sexcheck$IID %in% pedigree$MOTHER, "yes", NA)
sexcheck$Is.Dad <- ifelse(sexcheck$IID %in% pedigree$FATHER, "yes", NA)

#~~ cutoff between 0.6 and 0.8

sexcheck$NewSex <- ifelse(sexcheck$F < 0.6, "Female", ifelse(sexcheck$F > 0.8, "Male", "Unknown"))

table(sexcheck$NewSex, sexcheck$Is.Mum)
table(sexcheck$NewSex, sexcheck$Is.Dad)

#~~ Find bad genotyping IDs where Sex is dubious - these will be removed from the genabel dataset.

badgenoids <- c(sexcheck$id[which(sexcheck$NewSex == "Female" & sexcheck$Is.Dad == "yes")],
                sexcheck$id[which(sexcheck$NewSex == "Male" & sexcheck$Is.Mum == "yes")],
                sexcheck$id[which(sexcheck$NewSex == "Unknown")])

#~~ Now deal with sexes in newsheep

names(sexcheck)[2] <- "id"

#~~ extract family information from the PLINK file

system("cmd", input = paste0("gawk < ", plink.file, ".ped \"{print $1,$2,$3,$4,$5,$6}\" > plink.temp"))
famfile <- read.table("plink.temp")

names(famfile)[1:2] <- c("Family", "ID")

#~~ Add IDs that aren't in base file

length.marker <- length(famfile[!famfile$ID %in% newsheep$id,]$ID)

newsheep <- rbind(newsheep, data.frame(id      = famfile[!famfile$ID %in% newsheep$id,]$ID,
                                       NewHorn = rep(NA, length.marker),
                                       pedsex     = rep(0, length.marker),
                                       pedsex2    = rep("unknown", length.marker)))

newsheep <- join(newsheep, sexcheck[,c("id", "NewSex")])
head(newsheep)

newsheep$sex <- ifelse(newsheep$NewSex == "Male", 1, 0)
newsheep$sex[which(is.na(newsheep$sex))] <- 0

#~~ Create a horn phenotype vector for testing the data

newsheep$NewHorn2 <- NA
newsheep$NewHorn2[which(newsheep$NewHorn == 3 & newsheep$sex == 1)] <- 1
newsheep$NewHorn2[which(newsheep$NewHorn == 1 & newsheep$sex == 1)] <- 0

table(newsheep$id %in% famfile$ID)
table(famfile$ID %in% newsheep$id)

newsheep <- newsheep[newsheep$id %in% famfile$ID,]

#~~ create pheno file for GenABEL

write.table(newsheep, "data/GenABEL_hornpheno20150126.txt", sep = "\t", row.names = F, quote = F)

#~~ Create map file by extracting first 4 columns of tped file

mapfile <- read.table(paste0(plink.file,".map"))
names(mapfile)<- c("Chromosome", "SNP.Name", "MapDist", "Position")
mapfile <- mapfile[,-3]

#~~ Create GenABEL File

write.table(mapfile, paste0(plink.file, ".genabelmap"), col.names  = F, row.names = F, quote = F)

convert.snp.ped(pedfile = paste0(plink.file, ".ped"), 
                mapfile = paste0(plink.file, ".genabelmap"),
                outfile = "data/1_GenABEL_sheepabel_raw.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = F)   #added traits = 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Test GenABEL files                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sheepabel <- load.gwaa.data(phe = "data/GenABEL_hornpheno20150126.txt", 
                            gen = "data/1_GenABEL_sheepabel_raw.gen")

system("cmd", input = paste("plink --file", plink.file, "--sheep --freq"))

test.plink <- read.table("plink.frq", header = T)
head(test.plink)
names(test.plink)[3:4] <- c("A1.Plink", "A2.Plink")

test.plink <- subset(test.plink, CHR != 27)

test.genabel <- summary.snp.data(gtdata(sheepabel))
test.genabel$SNP <- row.names(test.genabel)
test.both <- join(test.plink, test.genabel)

plot(test.both$MAF, test.both$Q.2)

test.both[which((test.both$Q.2 - test.both$MAF) > 0.1),]

horntest <- qtscore(as.numeric(NewHorn2) ~ sex, sheepabel, trait = "binomial")

plot(horntest)
lambda(horntest)
summary(horntest)

# Top hit should be "OAR10_29511510.1"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Conduct quality control on GenABEL files            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove the badids 

nids(sheepabel)
if(length(badgenoids) > 0) sheepabel <- sheepabel[which(!idnames(sheepabel) %in% badgenoids),]

#~~ DEAL WITH X CHROMOSOME

sheepabel_x <- sheepabel[,which(chromosome(sheepabel) == 27)]

sheepabel <- recodeChromosome(sheepabel, list("27" = "X"), quiet = FALSE)
table(chromosome(sheepabel))


# We will use a perid.call of 0.99 for the time-being, and a minor allele frequency of 0.05, and an ibs.threshold 
# of 0.9. THe following is the output from R, showing the rounds of QC followed until all arguments are fulfilled:

qc.1to6 <- check.marker(sheepabel, perid.call = 0.99, maf = 0.01, ibs.threshold = 0.9)
save(qc.1to6, file = "results/qc_check_temp.Rdata")


#~~ use Xfix

sheepabel <- Xfix(data = sheepabel)



#~~ get and check QC'ed data

sheepabel2 <- sheepabel[qc.1to6$idok, qc.1to6$snpok]

#~~ Add the X pseudoautosomal region

qc.X <- check.marker(sheepabel_x, perid.call = 0.90, maf = 0.01, ibs.mrk = 0)
sheepabel_x <- sheepabel_x[,qc.X$snpok]

sheepabel_x <- recodeChromosome(sheepabel_x, list("27" = "X"), quiet = FALSE)

which(snp.names(sheepabel_x) %in% snp.names(sheepabel2))

nsnps(sheepabel_x)


#~~ Make sure the same IDs are coded in each dataset.

sheepabel2  <- sheepabel2 [intersect(idnames(sheepabel_x), idnames(sheepabel2)),]
sheepabel_x <- sheepabel_x[intersect(idnames(sheepabel_x), idnames(sheepabel2)),]

#~~ MERGE the datasets

sheepabel3 <- merge.gwaa.data(sheepabel_x, sheepabel2)
nsnps(sheepabel3)
nids(sheepabel3)
sheepabel3 <- sheepabel3[idnames(sheepabel2),]

if(length(badgenoids) > 0) sheepabel3 <- sheepabel3[which(!idnames(sheepabel3) %in% badgenoids),]


#~~ make sure everything corresponds with the old mapfile.

head(mapfile)

abelmap <- data.frame(SNP.Name = snp.names(sheepabel3),
                      Position.abel = map(sheepabel3),
                      Chr.abel = chromosome(sheepabel3))

testmap <- join(mapfile, abelmap, type = "inner")
head(testmap)

testmap$Chr.abel <- gsub("X", 27, testmap$Chr.abel)

which(testmap$Chromosome != testmap$Chr.abel)
testmap[which(testmap$Position   != testmap$Position.abel),]


#~~ Read in a file giving information on wrongly mapped SNPs

nsnps(sheepabel3)

#~~ Save dataset

sheepabel3 <- sheepabel3[,which(!snp.names(sheepabel3) %in% badsnps$SNP.Name)]
sheepabel3 <- sheepabel3[idnames(sheepabel3)[which(!is.na(phdata(sheepabel3)$NewSex))],]

save.gwaa.data(sheepabel3, phenofile = "data/1_GenABEL_hornpheno20150126.txt",
                           genofile  = "data/1_GenABEL_sheepabelQCed20150126.gen", human = FALSE)


#~~ if results directory doesn't exist, create it

if(!"results" %in% dir()) system("mkdir results")

write.table(qc.X$snpok, "results/1_Pseudoautosomal_SNPs_in_X.txt", row.names = F, quote = F)

system("rm plink.*")





# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # 5. Create a CRIMAP test dataset to deal with wrongly placed loci  #
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# # Take a subset of individuals as this is a very time consuming process.
# 
# #~~ remove non-genotyped parents from the pedigree
# 
# pedigree <- subset(pedigree, ANIMAL %in% idnames(sheepabel3))
# pedigree$MOTHER[which(!pedigree$MOTHER %in% idnames(sheepabel3))] <- 0
# pedigree$FATHER[which(!pedigree$FATHER %in% idnames(sheepabel3))] <- 0
# 
# #~~ remove parent if only one parent
# pedigree$MOTHER[which(pedigree$FATHER == 0 & pedigree$MOTHER != 0)] <- 0
# pedigree$FATHER[which(pedigree$FATHER != 0 & pedigree$MOTHER == 0)] <- 0
# 
# #~~ Find all offspring with two parents
# 
# offped <- pedigree[which(pedigree[,2] != 0),]
# offped$MumParents <- ifelse(offped$MOTHER %in% offped$ANIMAL, "yes", "no")
# offped$DadParents <- ifelse(offped$FATHER %in% offped$ANIMAL, "yes", "no")
# 
# offped <- subset(offped, MumParents == "yes" & DadParents == "yes")[,-c(4:5)]
# 
# #~~ create a Family Pedigree
# 
# famped <- NULL
# 
# for(i in offped$ANIMAL){
#   
#   ped1 <- offped[which(offped$ANIMAL == i),]
#   ped1 <- rbind(pedigree[which(pedigree$ANIMAL %in% c(ped1[,2:3])),], ped1)
#   
#   #~~ Make founder parents
#   
#   ped1 <- rbind(data.frame(ANIMAL = stack(ped1[1:2, 2:3])$values, FATHER = 0, MOTHER = 0), ped1)
#   ped1$Family <- paste("Offspring_", i, sep = "_")  
#   
#   famped <- rbind(famped, ped1)
# }
# 
# which(famped$FATHER %in% pedigree$MOTHER & famped$FATHER != 0)
# which(famped$MOTHER %in% pedigree$FATHER & famped$MOTHER != 0)
# 
# 
# #~~ Select the first 100 families
# 
# famped <- famped[which(famped$Family %in% unique(famped$Family)[1:100]),]
# 
# write.table(famped, 
#             "../results/1_FamPed_Raw_InitCrimapTest.txt", 
#             quote = F, sep = "\t", row.names = F)
# 
# #~~ Create Crimap Files for Autosomes
# 
# setwd("../crimap/crimap_test")
# 
# AnalysisSuffix <- "test"
# 
# system.time({
#   for(i in 1:26){
#     chrNo  <- genabelSubset(sheepabel3, i)
#     createCrimapFiles(chrNo, chromosomeid = paste(i, AnalysisSuffix, sep = ""), familyPedigree = famped)
#   }
# })
# 
# #~~ Create Crimap Files for X chromosome
# 
# chrNo  <- genabelSubset(sheepabel3, "X")
# createCrimapFiles(chrNo, chromosomeid = paste("X", AnalysisSuffix, sep = ""), familyPedigree=famped, pseudoautoSNPs = qc.X$snpok)
# system("cmd", input = paste0("rename chrX", AnalysisSuffix, ".gen chr27", AnalysisSuffix, ".gen"))
# 
# if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)
# 
# for(i in c(1:27)){
#   system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," prepare < crimapinput1 > chr", i, AnalysisSuffix, ".pre"))
# }
# 
# #~~ Run Chrompic and extract the linkage map
# 
# for(i in c(1:27)){
#   system("cmd", input = paste0("\"../crimap2504.exe\" ", i, AnalysisSuffix," chrompic > chr", i, AnalysisSuffix, ".cmp"))
#   system("cmd", input = "del *.cg")
# }
# 
# cmp.map <- NULL
# 
# for(i in 1:27) cmp.map <- rbind(cmp.map, MapFromChrompic(chrompicfile = paste0("chr", i, AnalysisSuffix, ".cmp")))
# 