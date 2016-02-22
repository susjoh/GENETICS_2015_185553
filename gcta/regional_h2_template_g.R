



recsumm <- read.table("../2_TotalSummIndivRR_FullCleanPostSim_g.txt", header = T)
rectab <- read.table("../2_IndivRR_FullCleanPostSim_g.txt", header = T)

parameters <- readLines("parameters.txt")

recsumm$RRID <- factor(recsumm$RRID)
recsumm$RRID2 <- recsumm$RRID

if(type == "cistrans"){
  recsumm$AnalysisRecombRate  <- recsumm$TotalRecombRate  
}

if(type == "trans"){
  
  rectab <- subset(rectab, Chr != chr)
  rectab$RecombCount.v2 <- as.numeric(rectab$RecombCount.v2)
  rectab$Inf.Chr.Length <- as.numeric(rectab$Inf.Chr.Length)
  
  tempsumm <- data.frame(AnalysisRecombCount = unlist(tapply(rectab$RecombCount.v2, rectab$UniqueID2, sum)),
                         AnalysisInfLength = unlist(tapply(rectab$Inf.Chr.Length, rectab$UniqueID2, sum)))
  tempsumm$UniqueID2 <- row.names(tempsumm)
  tempsumm$AnalysisRecombRate <- (tempsumm$AnalysisRecombCount/tempsumm$AnalysisInfLength) * 1e6
  recsumm <- merge(recsumm, tempsumm, by = "UniqueID2")
}

if(type == "cis"){
  
  rectab <- subset(rectab, Chr == chr)
  rectab$RecombCount.v2 <- as.numeric(rectab$RecombCount.v2)
  rectab$Inf.Chr.Length <- as.numeric(rectab$Inf.Chr.Length)
  
  tempsumm <- data.frame(AnalysisRecombCount = rectab$RecombCount.v2,
                         AnalysisInfLength = rectab$Inf.Chr.Length,
                         UniqueID2 = rectab$UniqueID2)
  tempsumm$AnalysisRecombRate <- (tempsumm$AnalysisRecombCount/tempsumm$AnalysisInfLength) * 1e6
  recsumm <- merge(recsumm, tempsumm, by = "UniqueID2")
}

recsumm.m <- droplevels(subset(recsumm, RRID.SEX == "Male"))
recsumm.f <- droplevels(subset(recsumm, RRID.SEX == "Female"))

library(asreml)

source("../makeGRM.R")
source("../ASReml.EstEffects.R")
source("../ASReml.ExtractPredictors.R")


# ALL IDS

grm.rest <- read.table(paste0("150204_wo",i,"_GRM_adj.grm.gz"))
ids.rest <- read.table(paste0("150204_wo",i,"_GRM_adj.grm.id"))

grm.reg <- read.table(paste0("150204_",i,"_GRM_adj.grm.gz"))
ids.reg <- read.table(paste0("150204_",i,"_GRM_adj.grm.id"))

reginv  <- makeGRM(grm.reg, ids.reg)
restinv <- makeGRM(grm.rest, ids.rest)


model1 <- asreml(fixed = AnalysisRecombRate ~ factor(RRID.SEX) + RRID.Fhat3,
                 random = ~ giv(RRID) + ide(RRID),
                 data = recsumm,
                 ginverse =  list(RRID = restinv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6)

model2 <- asreml(fixed = AnalysisRecombRate ~ factor(RRID.SEX) + RRID.Fhat3,
                 random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                 data = recsumm,
                 ginverse =  list(RRID = restinv, RRID2 = reginv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6)

obj1.1 <- summary(model1, all = T)$coef.fixed
obj1.2 <- ASReml.EstEffects(model1)
obj1.3 <- summary(model1, all = T)$loglik 

obj2.1 <- summary(model2, all = T)$coef.fixed
obj2.2 <- ASReml.EstEffects(model2)
obj2.3 <- summary(model2, all = T)$loglik 

obj1.1
obj1.2
obj1.3

obj2.1
obj2.2
obj2.3

rm(model1, model2)
gc()

model1.m <- asreml(fixed = AnalysisRecombRate ~ RRID.Fhat3,
                   random = ~ giv(RRID) + ide(RRID),
                   data = recsumm.m,
                   ginverse =  list(RRID = restinv),
                   na.method.X = "omit", na.omit.Y = "na.omit",
                   workspace = 500e+6, pworkspace = 500e+6)

model2.m <- asreml(fixed = AnalysisRecombRate ~ RRID.Fhat3,
                   random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                   data = recsumm.m,
                   ginverse =  list(RRID = restinv, RRID2 = reginv),
                   na.method.X = "omit", na.omit.Y = "na.omit",
                   workspace = 500e+6, pworkspace = 500e+6)

obj1.1.m <- summary(model1.m, all = T)$coef.fixed
obj1.2.m <- ASReml.EstEffects(model1.m)
obj1.3.m <- summary(model1.m, all = T)$loglik 

obj2.1.m <- summary(model2.m, all = T)$coef.fixed
obj2.2.m <- ASReml.EstEffects(model2.m)
obj2.3.m <- summary(model2.m, all = T)$loglik 

obj1.1.m
obj1.2.m
obj1.3.m

obj2.1.m
obj2.2.m
obj2.3.m

rm(model1.m, model2.m)
gc()

model1.f <- asreml(fixed = AnalysisRecombRate ~ RRID.Fhat3,
                   random = ~ giv(RRID) + ide(RRID),
                   data = recsumm.f,
                   ginverse =  list(RRID = restinv),
                   na.method.X = "omit", na.omit.Y = "na.omit",
                   workspace = 500e+6, pworkspace = 500e+6)

model2.f <- asreml(fixed = AnalysisRecombRate ~ RRID.Fhat3,
                   random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                   data = recsumm.f,
                   ginverse =  list(RRID = restinv, RRID2 = reginv),
                   na.method.X = "omit", na.omit.Y = "na.omit",
                   workspace = 500e+6, pworkspace = 500e+6)



obj1.1.f <- summary(model1.f, all = T)$coef.fixed
obj1.2.f <- ASReml.EstEffects(model1.f)
obj1.3.f <- summary(model1.f, all = T)$loglik 

obj2.1.f <- summary(model2.f, all = T)$coef.fixed
obj2.2.f <- ASReml.EstEffects(model2.f)
obj2.3.f <- summary(model2.f, all = T)$loglik 

obj1.1.f
obj1.2.f
obj1.3.f

obj2.1.f
obj2.2.f
obj2.3.f

rm(model1.f, model2.f)
gc()


save(obj1.1  , obj1.2  , obj1.3  , obj2.1  , obj2.2  , obj2.3  ,
     obj1.1.m, obj1.2.m, obj1.3.m, obj2.1.m, obj2.2.m, obj2.3.m,
     obj1.1.f, obj1.2.f, obj1.3.f, obj2.1.f, obj2.2.f, obj2.3.f,file = paste0("model", i, ".RData"))
