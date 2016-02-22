# Functions for calculating recombination counts
# Author: Susan Johnston


#~~ First define functions to determine recombination Count of the chromosome

RecCountFunc <- function(vector){
  
  unlist(lapply(vector, function(test){
    
    y <- unlist(strsplit(test, split = ""))
    y <- y[which(y != "-")]
    y[which(y %in% c("o", "c"))] <- 0
    y[which(y %in% c("i", ":"))] <- 1
    
    y1 <- length(which(c(-9, y) != c(y, -9))) - 2
    if(y1 == -2) y1 <- 0
    
    return(y1)
  }
  )
  )
}

ConservativeRecCountFunc <- function(vector){
  
  unlist(lapply(vector, function(test){
    
    y <- unlist(strsplit(test, split = ""))
    y <- y[which(!y %in% c("-", ":", "c"))]
    y[which(y == "o")] <- 0
    y[which(y == "i")] <- 1
    
    return(length(which(c(-9, y) != c(y, -9))) - 2)
  }
  )
  )
}



# Here we determine which SNP came first in ORDER (i.e. SNP 1, SNP 2, etc.~)


InfLengthFunc <- function(CMPstring, position = "First"){
  test <- gsub(" ", "", CMPstring)
  temp <- which(unlist(strsplit(test, split = "")) %in% c("i", "1","0","o", ":", "c"))
  if(position == "First") return(temp[1])
  if(position == "Last"){
    if(length(temp) == 0) return(NA)
    if(length(temp) > 0) return(temp[length(temp)])
  }
}  


InfLengthFunc.Cons <- function(CMPstring, position = "First"){
  test <- gsub(" ", "", CMPstring)
  temp <- which(unlist(strsplit(test, split = "")) %in% c("i", "1","0","o"))
  if(position == "First") return(temp[1])
  if(position == "Last"){
    if(length(temp) == 0) return(NA)
    if(length(temp) > 0) return(temp[length(temp)])
  }
}  