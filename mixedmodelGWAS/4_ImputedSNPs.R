

library(GenABEL)
library(asreml)
library(plyr)

impabel <- load.gwaa.data(genofile = "imputedPLINK.genabel", phenofile = "imputedPLINK.phenofile")


#chr6results <- read.table("../mach/ImputationDataForPlotting.txt", header = T)



load("GLMM.data.Rdata")


wald.results <- NULL
effect.results <- NULL

analysis.type <- "cistrans"  #"cistrans"

rectab <- droplevels(subset(rectab, UniqueID2 %in% recsumm$UniqueID2))

it.start <- 1
it.end <- nsnps(impabel)

for(i in it.start:it.end){
  
  if(i %in% seq(1, it.end, 10)) print(paste("Analysing SNP", i, "of", it.end))
  
  if(analysis.type == "cistrans"){
    temprecsumm <- data.frame(as.character.gwaa.data(sheepabel[,i]))
    temprecsumm$RRID <- row.names(temprecsumm)
    names(temprecsumm)[1] <- "SNP"
    temprecsumm <- join(newrecsumm, temprecsumm)
  }
  
  if(analysis.type == "trans"){
    
    temprectab <- subset(rectab, Chr != chromosome(sheepabel)[i])
    temprecsumm <- data.frame(TotalRecombCount = tapply(temprectab$RecombCount.v2, temprectab$UniqueID2, sum))
    temprecsumm$RRID <- sapply(row.names(temprecsumm), function(x) gsub("RRID", "", strsplit(x, split = "_")[[1]][4]))
    
    
    temprecsumm1 <- data.frame(as.character.gwaa.data(sheepabel[,i]))
    temprecsumm1$RRID <- row.names(temprecsumm1)
    names(temprecsumm1)[1] <- "SNP"
    
    temprecsumm <- join(temprecsumm, temprecsumm1)
    temprecsumm <- join(temprecsumm, unique(subset(newrecsumm, select = -TotalRecombCount)))
    head(temprecsumm)
    rm(temprecsumm1, temprectab)
    temprecsumm$RRID <- as.factor(temprecsumm$RRID)
    
  }
  
  
  fitall <- asreml(fixed = TotalRecombCount ~ RRID.SEX + RRID.Fhat3 + as.factor(SNP),
                   random = ~ giv(RRID) + ide(RRID),
                   data = temprecsumm,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit",
                   workspace = 500e+6, pworkspace = 500e+6, trace = F)
  
  fit.m <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3 + as.factor(SNP),
                  random = ~ giv(RRID) + ide(RRID),
                  data = droplevels(subset(temprecsumm, RRID.SEX == "Male")),
                  ginverse =  list(RRID = ainv),
                  na.method.X = "omit", na.omit.Y = "na.omit",
                  workspace = 500e+6, pworkspace = 500e+6, trace = F)
  
  fit.f <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3 + as.factor(SNP),
                  random = ~ giv(RRID) + ide(RRID),
                  data = droplevels(subset(temprecsumm, RRID.SEX == "Female")),
                  ginverse =  list(RRID = ainv),
                  na.method.X = "omit", na.omit.Y = "na.omit",
                  workspace = 500e+6, pworkspace = 500e+6, trace = F)
  
  
  
  x1 <- data.frame(summary(fitall, all = T)$coef.fixed)
  x2 <- data.frame(summary(fit.m , all = T)$coef.fixed)
  x3 <- data.frame(summary(fit.f , all = T)$coef.fixed)
  
  x1$Effect <- row.names(x1)
  x2$Effect <- row.names(x2)
  x3$Effect <- row.names(x3)
  
  x4 <- rbind(cbind(x1, Model = "All"),
              cbind(x2, Model = "Male"),
              cbind(x3, Model = "Female"))
  
  x4 <- x4[-grep("RRID", x4$Effect),]
  x4$SNP.Name <- snpnames(sheepabel)[i]
  
  rm(x1, x2, x3)
  
  x1 <- data.frame(wald.asreml(fitall))
  x2 <- data.frame(wald.asreml(fit.m))
  x3 <- data.frame(wald.asreml(fit.f))
  x1 <- x1[grep("SNP", row.names(x1)),]
  x2 <- x2[grep("SNP", row.names(x2)),]
  x3 <- x3[grep("SNP", row.names(x3)),]
  
  x5 <- rbind(cbind(x1, Model = "All"),
              cbind(x2, Model = "Male"),
              cbind(x3, Model = "Female"))
  
  x5$SNP.Name <- snpnames(sheepabel)[i]
  
  wald.results <- rbind(wald.results, x5)
  effect.results <- rbind(effect.results, x4)
  rm(x1, x2, x3, x4, x5, fitall, fit.m, fit.f, temprecsumm)
  
}

write.table(wald.results, paste0("ImputedWaldtest.It_", it.start, "_", it.end, "_", analysis.type, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(effect.results, paste0("ImputedEffectRes.It_", it.start, "_", it.end, "_", analysis.type, ".txt"), row.names = F, sep = "\t", quote = F)
