

no.sim <- 100


for(i in 1:no.sim){
  
  system(paste0("mkdir sim", i))
  
  writeLines(paste0("AnalysisSuffix <- \"sim", i, "\""), paste0("sim", i, "/simscript.R"))
  system(paste0("cat 1_Simulate_Recombination.R >> sim", i, "/simscript.R"))
  
  
  writeLines(paste0("#!/bin/sh\n", 
	"\n",
	"#$ -cwd\n",
	"#$ -l h_rt=12:00:00\n",
	"#$ -V\n",
	"#$ -l h_vmem=5200M\n",
	"\n",
	". /etc/profile.d/modules.sh\n",
	"module load R/3.1.1\n",
	"R CMD BATCH simscript.R\n"),
	paste0("sim", i, "/sim.sh"))
             

             
  
  setwd(paste0("sim", i))
  system(paste0("qsub sim.sh"))
  setwd("..")
  
}
