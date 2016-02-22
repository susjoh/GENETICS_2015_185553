#~~ Define function

for(i in 1:1){
  
  write(paste("iteration = ", i), file = paste0("ExtractRecsumm", i, ".R"))
  
  system(paste0("cat 2_AnalysisSimulationOutput.R >> ExtractRecsumm", i, ".R"))
  
    writeLines(paste0("#!/bin/sh\n", 
                      "\n",
                      "#$ -cwd\n",
                      "#$ -l h_rt=47:00:00\n",
                      "#$ -V\n",
                      "#$ -l h_vmem=5200M\n",
                      "\n",
                      ". /etc/profile.d/modules.sh\n",
                      "module load R/3.1.1\n",
                      "R CMD BATCH ExtractRecsumm", i, ".R\n"),
               paste0("ExtractRecsumm", i, ".sh"))
    
    system(paste0("qsub ExtractRecsumm", i, ".sh"))
}
    




