#
# GCTA Analysis
# Susan Johnston
# June 2014
#

library(plyr)

# Analyses outlines

# prefix a: 150 SNP windows overlapping trans
# prefix b: 50 SNP windows overlapping trans
# prefix 

#~~  Specify Analysis Parameters

prefix <- "d.ps"
sw <- 20	
overlap <- 10
chrvec <- 7	
chrvecname <- "7"
type   <- "cistrans"         #trans, cis, cistrans
markersubset <- NULL
positionsubset <- NULL #c(110000000, 118000000)
#sex <- "both"    #"M" or "F" or "both"
plink.prefix <- "20150129merged1_66nodups.QC2"

#~~ Create working directory

analysis_name <- paste0("Run_", prefix, "_", sw, "SNP_Chr", chrvecname, "_", type)

snplistvec <- NULL

if(!analysis_name %in% dir()) system(paste("mkdir", analysis_name))

setwd(analysis_name)


#~~ Extract the first 2 columns from the .ped file

# system(paste0("cut -d\'*\' -f1-2  ../", plink.prefix, ".ped > idlist.txt"))




#~~ Generate the SNP list files

mapfile <- read.table(paste0("../", plink.prefix, ".map"), stringsAsFactors = F)
head(mapfile)

names(mapfile) <- c("Chr", "SNP.Name", "cM", "Position")
mapfile <- arrange(mapfile, Chr, Position)


for(i in chrvec){
  
  x <- subset(mapfile, Chr == i)$SNP.Name
  
  if(!is.null(markersubset)) x <- x[markersubset,]
  if(!is.null(positionsubset)){
  	x <- subset(mapfile, Chr == i & Position > positionsubset[1] & Position < positionsubset[2])$SNP.Name
   }
  
  for(j in seq(1, length(x), overlap)){
    
    iteration_code <- paste0("chr", i, ".", prefix, ".", j,".",sw)
    
    write.table(na.omit(x[j:(j + (sw-1))]),
                paste0("chr", i, ".", prefix, ".", j,".",sw, "snplist.txt"),
                row.names = F, col.names = F, quote = F)
    
    snplistvec <- c(snplistvec, iteration_code)
    
  }
  
  
  iteration_code <- paste0("chr", i, ".", prefix, ".", (length(x)-sw+1),".",sw)
  
  write.table(na.omit(x[(length(x)-sw+1):length(x)]),
              paste0("chr", i, ".", prefix, ".", (length(x)-sw+1),".",sw, "snplist.txt"),
              row.names = F, col.names = F, quote = F)
  
  snplistvec[length(snplistvec)]  <- iteration_code
  
}

write.table(snplistvec, "snplistvec.txt", row.names = F, col.names = F, quote = F)

write.table(c(prefix, sw, overlap, paste(chrvec, collapse = " "), chrvecname, type, analysis_name),
            "parameters.txt", row.names = F, col.names = F, quote = F)


system("cp ../Regional_heritability_asreml.R Regional_heritability_asreml.R")

system("R CMD BATCH Regional_heritability_asreml.R")

