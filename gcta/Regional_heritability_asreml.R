
snplistvec <- read.table("snplistvec.txt", colClasses = "character")
parameters <- readLines("parameters.txt")

### Run each individual chromosome

for(i in snplistvec$V1){

  chr <- gsub("chr", "", strsplit(i, "\\.")[[1]][1])
  type <- parameters[6]
  
  write.table(paste0("i <- \"",i, "\""), paste0("partition_", i, ".R"), quote = F, row.names = F, col.names = F)
  write.table(paste0("chr <- \"",chr, "\""), paste0("partition_", i, ".R"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("type <- \"",type, "\""), paste0("partition_", i, ".R"), quote = F, row.names = F, col.names = F, append = T)
  
  
  system(paste0("cat ../regional_h2_template_g.R partition_", i, ".R >> partition_", i, ".R"))
  
  
  write.table(paste0("#!/bin/sh"         ), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F)
  write.table(paste0(""                  ), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("#$ -cwd"           ), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("#$ -l h_rt=00:30:00"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("#$ -V"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("#$ -l h_vmem=5200M"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0(""                  ), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0(". /etc/profile.d/modules.sh"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("module load R/2.15.2"                  ), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  
  
  write.table(paste0("../gcta64 --bfile ../20150129merged1_66nodups.QC2 --autosome-num 26 --extract ", i, "snplist.txt --keep ../idlist.txt --make-grm-gz --out 150204_", i, "_GRM"),
              paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("../gcta64 --grm-gz 150204_", i, "_GRM --grm-adj 0 --make-grm-gz --out 150204_", i, "_GRM_adj"),
              paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  
  
  write.table(paste0("../gcta64 --bfile ../20150129merged1_66nodups.QC2 --autosome-num 26 --autosome --exclude ", i, "snplist.txt --keep ../idlist.txt --make-grm-gz --out 150204_wo", i, "_GRM"),
	      paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste0("../gcta64 --grm-gz 150204_wo", i, "_GRM --grm-adj 0 --make-grm-gz --out 150204_wo", i, "_GRM_adj"),
              paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
  
  write.table(paste0("R CMD BATCH partition_", i, ".R"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)

#  write.table(paste0("rm 150204_", i, "_*"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)
#  write.table(paste0("rm 150204_wo", i, "_*"), paste0("reg", i ,".sh"), quote = F, row.names = F, col.names = F, append = T)

  if(i != nrow(snplistvec)) system(paste0("qsub reg", i ,".sh"))
  if(i == nrow(snplistvec)) system(paste0("qsub reg", i ,".sh -m e -M sjohns10@staffmail.ed.ac.uk"))

}

