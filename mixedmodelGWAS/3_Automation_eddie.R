# 
start.no <- 1
stop.no <- 39418
analysis.count <- 100

template <- readLines("2_MMGWAS_Template.R")

for(i in seq(start.no, stop.no, analysis.count)){
  stopval <- (i + analysis.count - 1)
  if(stopval > stop.no) stopval <- stop.no
  
  writeLines(paste("it.start <- ", i, "\n",
                   "it.end <- ", stopval, "\n"),
             paste0("it", i, "_", stopval, ".R"))
  write.table(template, paste0("it", i, "_", stopval, ".R"), append = T, row.names = F, quote = F, col.names = F)
  
  writeLines(paste0("#!/bin/sh\n", 
                    "\n",
                    "#$ -cwd\n",
                    "#$ -l h_rt=01:00:00\n",
                    "#$ -V\n",
                    "#$ -l h_vmem=5200M\n",
                    "\n",
                    ". /etc/profile.d/modules.sh\n",
                    "module load R/3.1.1\n",
                    "R CMD BATCH it", i, "_", stopval, ".R\n"),
             paste0("it", i, "_", stopval, ".sh"))
  
  system(paste0(" qsub it", i, "_", stopval, ".sh"))
}

