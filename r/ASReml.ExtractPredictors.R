#~~ Extract predictors from ASReml models


extractPredictors <- function(model, idvec, idvar){
  
  #~~ Extract the residuals and Linear Predictors and calculate the mean per individual
  
  residtab <- data.frame(ID       = as.character(idvec),
                         Residuals  = as.numeric(model$residuals),
                         LinearPred = as.numeric(model$linear.predictors))
  
  residtab <- data.frame(Resid.Mean = tapply(residtab$Residuals,  residtab$ID, mean, na.rm = T),
                         LPred.Mean = tapply(residtab$LinearPred, residtab$ID, mean, na.rm = T))
  
  residtab$ID <- row.names(residtab)
  
  attributes(residtab$Resid.Mean) <- NULL
  attributes(residtab$LPred.Mean) <- NULL
  
  
  #~~ Extract BLUPs for the individual
  
  bluptab        <- data.frame(BLUP = model$coefficients$random)
  bluptab$ID     <- row.names(bluptab)
  bluptab$Effect <- unlist(lapply(bluptab$ID, function (x) strsplit(x, split = "_")[[1]][1]))
  bluptab$ID     <- unlist(lapply(bluptab$ID, function (x) strsplit(x, split = "_")[[1]][2]))
  
  bluptab <- bluptab[grep(idvar, bluptab$Effect),]
  bluptab <- cast(bluptab, value= "BLUP", ID ~ Effect)
  
  names(bluptab)[2:ncol(bluptab)] <- paste(substr(names(bluptab)[2:ncol(bluptab)], 1, 3), "BLUP", sep = ".")
  
  #~~ join into a single results table
  
  require(plyr)
  
  predtab <- join(bluptab, residtab, type = "right")
  
  names(predtab)[which(names(predtab) == "ID")] <- idvar
  
  predtab
  
}