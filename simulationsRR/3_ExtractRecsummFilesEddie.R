for(i in c(1:100)){
  system(paste0("cp sim", i, "/SummaryModels", i, ".Rdata recsummfiles/SummaryModels", i, ".Rdata"))
}

  