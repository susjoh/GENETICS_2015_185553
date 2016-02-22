# Functions for extracting information from Crimap Chrompic Files
# Author: Susan Johnston

#~~ GGPLOT2 default options

library(ggplot2)

defopts <- theme(axis.text.x  = element_text (size = 16, vjust = 0),
                 axis.text.y  = element_text (size = 14, hjust = 1.3),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 strip.background = element_blank())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Pull Map from Chrompic                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


MapFromChrompic <- function(chrompicfile){
  
  x <- readLines(chrompicfile)
  map <- x[grep("Sex_averaged", x):length(x)]
  map <- map[3:length(map)]
  map <- map[1:(length(map)-6)]
  
  map <- data.frame(map, stringsAsFactors=F)
  
  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]  
  
  for(i in 1:10) map$map  <- gsub("  ", " ", map$map)
  for(i in 1:10) map$map2 <- gsub("  ", " ", map$map2)
  
  rm(x)
  
  map$Order <- NA
  map$SNP.Name <- NA
  map$Position <- NA
  map$r <- NA
  map$cMdiff <- NA
  
  for(i in 1:nrow(map)){
    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]
    
    map$Order[i] <- x[1]
    map$SNP.Name[i] <- x[2]
    map$Position[i] <- x[3]
    rm(x)
    
    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]
    
    map$r[i] <- x[1]
    map$cMdiff[i] <- x[2]
    
    rm(x)
  }
  
  map <- subset(map, select = -c(map, map2))
  
  map
  
  
}



MapFromMap <- function(mapfile){
  
  x <- readLines(mapfile)
  map <- x[grep("Sex-specific map", x):length(x)]
  
  
  map <- map[3:length(map)]
  map <- map[1:(length(map)-7)]
  map <- map[-which(map == "")]
  
  map <- data.frame(map, stringsAsFactors=F)
  
  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]  
  
  for(i in 1:1e4) if(length(grep("  ", map$map)  > 0)) map$map   <- gsub("  ", " ", map$map)  else break
  for(i in 1:1e4) if(length(grep("  ", map$map2) > 0)) map$map2  <- gsub("  ", " ", map$map2) else break
  
  rm(x)
  
  map$Order <- NA
  map$SNP.Name <- NA
  map$cMPosition.Female <- NA
  map$cMPosition.Male <- NA
  map$Female.r <- NA
  map$cMdiff.Female <- NA
  map$Male.r <- NA
  map$cMdiff.Male <- NA
  
  for(i in 1:nrow(map)){
    
    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]
    
    map$Order[i]             <- x[1]
    map$SNP.Name[i]          <- x[2]
    map$cMPosition.Female[i] <- x[3]
    map$cMPosition.Male[i]   <- x[4]
    
    rm(x)
    
    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]
    
    map$Female.r[i]      <- x[1]
    map$cMdiff.Female[i] <- x[2]
    map$Male.r[i]        <- x[3]
    map$cMdiff.Male[i]   <- x[4]
    
    
    rm(x)
  }
  
  map <- subset(map, select = -c(map, map2))
  
  map
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Determine the recombination rate form chrompic files       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


RecombRateFromChrompic <- function(chrompicfile){
  
  #~~ read lines from the chrompic file
  
  x <- readLines(chrompicfile)
  
  nmarkers <- nrow(MapFromChrompic(chrompicfile))
  
  x <- x[1:grep("Sex_averaged", x)]
  
  rowstokeep <- sort(c(grep("^ ", x),
                       grep("phase likelihood", x)))
  
  x <- x[rowstokeep]
  if(nmarkers < 101) x <- x[-c(1:nmarkers)]
  if(nmarkers > 100) x <- x[-c(1:100)]
  
  
  x <- x[-grep("         ", x)]
  
  #~~ create a table with family information to add to the full table later
  
  famtab <- data.frame(FamilyOrder = 1:length(grep("Family", x)),
                       FamilyID = grep("Family", x, value = T),
                       LineNo = grep("Family", x),
                       stringsAsFactors = F)
  
  famtab$startline <- famtab$LineNo - (famtab$FamilyOrder - 1)
  famtab$FamilyShortID <- NA
  for(i in 1:nrow(famtab)) famtab$FamilyShortID[i] <- strsplit(famtab$FamilyID[i], split = " ")[[1]][2]
  
  x <- x[-grep("Family", x)]
  
  famtab$Count <- diff(c(famtab$startline, (length(x) + 1)))
  
  famvec <- rep(famtab$FamilyShortID, famtab$Count)
  
  #~~ create a data frame to put all information in
  
  recombframe <- data.frame(data = x, 
                            id = NA, 
                            RecombCount = NA,
                            parent = NA,
                            Family = famvec,
                            stringsAsFactors=F)
  
  
  recombframe$data <- gsub("^  ", " ",  recombframe$data)
  recombframe$data <- gsub("^   ", " ",  recombframe$data)
    
  #~~ Get the IDs
  
  
  recombframe$id <- unlist(lapply(recombframe$data, function(foo) strsplit(foo, split = " ")[[1]][2]))
  recombframe$parent[which(recombframe$id != "")] <- "MOTHER"
  recombframe$id[which(recombframe$id == "")] <- recombframe$id[which(recombframe$id != "")]
  recombframe$parent[which(is.na(recombframe$parent))] <- "FATHER"
  
  
  
  #~~ get the number of recombinations for each individual
  
  recExtract <- function(CMPstring){
    x <- strsplit(CMPstring, split = " ")[[1]][length(strsplit(CMPstring, split = " ")[[1]])]
  }

  recombframe$RecombCount <- sapply(recombframe$data, recExtract, simplify = TRUE)
    
  #~~ format the data column so that we can get out the number of informative loci
  removeID <- function(CMPstring){
    
    x <- strsplit(CMPstring, split = " ")[[1]]
    paste0(x[3:(length(x)-1)], collapse = "")

  }  

  system.time(recombframe$data <- unlist(lapply(recombframe$data, removeID)))

  recombframe$data2 <- gsub("-" , "", recombframe$data)
  recombframe$data2 <- gsub("c" , "", recombframe$data2)
  recombframe$data2 <- gsub(":" , "", recombframe$data2)
  
  recombframe$No.Inf.Loci <- nchar(recombframe$data2)
  
  recombframe <- subset(recombframe, select = -data2)
  
  head(recombframe)
  
  recombframe$RecombCount <- as.numeric(recombframe$RecombCount)
  
  recombframe
  
}



