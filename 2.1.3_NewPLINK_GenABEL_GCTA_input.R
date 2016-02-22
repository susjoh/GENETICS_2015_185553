# Prepare PLINK files for analysis in GCTA
# Author: Susan Johnston

# Requires plink, gawk, mkdir in PATH

library(plyr)
library(GenABEL)

#~~ Load sheepabel

sheepabel <- load.gwaa.data(phenofile = "data/1_GenABEL_hornpheno20150126.txt",
                            genofile =  "data/1_GenABEL_sheepabelQCed20150126.gen")

#~~ Load PAR snps

PAR.snps <- readLines("results/1_Pseudoautosomal_SNPs_in_X.txt")

#~~ Create a PLINK file with all IDs and SNPs passing QC in 1.0_Create_Genabel_Initial_QC.R

#~~ write SNP IDs to file
write.table(snp.names(sheepabel), "data/2_Final_SNPs.txt", row.names = F, quote = F, col.names = F)

#~~ Write IDs to file

system("cmd", input = paste0("gawk < data/20130410merged1_66nodups.oar3.1.ped \"{print $1,$2}\" > data/2_Final_IDs.txt"))
id.file <- read.table("data/2_Final_IDs.txt")
head(id.file)

id.file <- subset(id.file, V2 %in% idnames(sheepabel))
write.table(id.file,   "data/2_Final_IDs.txt" , row.names = F, quote = F, col.names = F)

system("plink --file data\\20130410merged1_66nodups.oar3.1 --sheep --recode --make-bed --keep data\\2_Final_IDs.txt --extract data\\2_Final_SNPs.txt --out data\\20150129merged1_66nodups.QC1")

#~~ revise to remove all Mendelian errors and heterozygous SNPs

menderr <- read.table("results/2_MendelianErrors_Crimap_g.txt", header = T, stringsAsFactors = F)
head(menderr)
menderr$Order <- menderr$Locus + 1

snpmap <- read.table("results/2_Maps_from_Chrompic_g.txt", header = T, stringsAsFactors = F)
head(snpmap)

menderr <- unique(join(menderr, snpmap[,c("Order", "SNP.Name", "Chr")]))
table(menderr$Chr)

famids <- read.table("data/2_Final_IDs.txt")
names(famids) <- c("FAMILY", "ANIMAL")

menderr <- join(menderr, famids)
head(menderr)

menderr <- menderr[,c("FAMILY", "ANIMAL", "SNP.Name")]
menderr <- cbind("Q", menderr, 0.01)
head(menderr)

names(menderr) <- c("Q", "FID", "IID", "SNPID", "score")

#~~ now find the heterozygous males from the previous PLINK analysis

hetmales <- read.table("plink.hh", stringsAsFactors = F)
head(hetmales)

hetmales <- subset(hetmales, !V3 %in% PAR.snps)
head(hetmales)

#~~ remove those that are actually female
hetmales <- hetmales[-which(hetmales$V2 %in% phdata(sheepabel)$id[which(phdata(sheepabel)$sex == 0)]),]
head(hetmales)

#~~ remove those from loci that were removed because of quality control

hetmales <- subset(hetmales, V3 %in% snp.names(sheepabel))

hetmales <- cbind("Q", hetmales, 0.01)
names(hetmales) <- c("Q", "FID", "IID", "SNPID", "score")

menderr <- rbind(menderr, hetmales)

write.table(menderr, "data/2_MendErr_Qscores.txt", row.names = F, col.names = F, sep = "\t", quote = F)

#~~ write an updated sex file

sextab <- phdata(sheepabel)[,c("id", "sex")]
names(sextab) <- c("ANIMAL", "SEX")

sextab <- join(sextab, famids)
head(sextab)

sextab$SEX[which(sextab$SEX == 0)] <- 2

sextab <- sextab[,c("FAMILY", "ANIMAL", "SEX")]

write.table(sextab, "data/2_NewSexForPlink.txt", row.names = F, col.names = F, sep = "\t", quote = F) 


#~~ RUN PLINK AGAIN

system("plink --bfile data/20150129merged1_66nodups.QC1 --make-bed --sheep --recode --qual-geno-scores data/2_MendErr_Qscores.txt --qual-geno-threshold 0.2 --update-sex data/2_NewSexForPlink.txt --out data/20150129merged1_66nodups.QC2")



table(hetmales$V2)

#~~ Output SNP lists for GCTA

if(!"gcta" %in% dir()) system("cmd", input = "mkdir gcta")

for(i in 0:26){
  write.table(snp.names(sheepabel)[chromosome(sheepabel) == i],
              paste0("gcta/chr", i, "snplist.txt"),
              row.names = F, col.names = F, quote = F)
}

table(chromosome(sheepabel))

write.table(snp.names(sheepabel)[chromosome(sheepabel) == "X"],
            "gcta/chr27snplist.txt",
            row.names = F, col.names = F, quote = F)

write.table(snp.names(sheepabel)[chromosome(sheepabel) == "X"][which(!snp.names(sheepabel)[chromosome(sheepabel) == "X"] %in% PAR.snps)],
            "gcta/chr27nonparsnplist.txt",
            row.names = F, col.names = F, quote = F)

write.table(PAR.snps,
            "gcta/chr27parsnplist.txt",
            row.names = F, col.names = F, quote = F)

  
#~~ Create a new GenABEL input with the cleaned SNP data

system("cmd", input = "plink --bfile data/20150129merged1_66nodups.QC2 --sheep --recode --out data/20150129merged1_66nodups.QC2")

mapfile <- read.table("data/20150129merged1_66nodups.QC2.map")
head(mapfile)

write.table(mapfile[,c(1, 2, 4)], "data/20150129merged1_66nodups.QC2.genabelmap",
            row.names = F, col.names = F, sep = "\t", quote = F)


convert.snp.ped(pedfile = "data/20150129merged1_66nodups.QC2.ped", 
                mapfile = "data/20150129merged1_66nodups.QC2.genabelmap",
                outfile = "data/1_GenABEL_sheepabelFullQC20150129.gen",
                strand = "u", bcast = 1000000, traits = 1, mapHasHeaderLine = F)   #added traits = 1



sheepabel <- load.gwaa.data(phenofile = "data/1_GenABEL_hornpheno20150126.txt",
                            genofile =  "data/1_GenABEL_sheepabelFullQC20150129.gen")

summary(qtscore(NewHorn2, data = sheepabel))
table(chromosome(sheepabel))
plot(map(sheepabel)[chromosome(sheepabel) == 27])

string <- readLines("data/20150129merged1_66nodups.QC2.ped", n = 1)
stringvec <- strsplit(string, split = " ")[[1]]
length(stringvec)/2

