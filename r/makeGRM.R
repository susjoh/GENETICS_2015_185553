# # Provided by Camillo Berenos
# 
# # makeGRM <- function(grm.object, ids.object, phenoframe, idcol = 1) {
# #   
# #   names(phenoframe)[idcol] <- "ANIMALx"
# #   phenoframe$ANIMALx <- as.factor(phenoframe$ANIMALx)
# #   
# #   elements <- length(ids.object[, 1])
# #   X <- diag(elements)
# #   
# #   X[upper.tri(X, diag = TRUE)] <- grm.object[, 4]
# #   X.grm <- X + t(X) - diag(diag(X))
# #   nam <- as.factor(ids.object[, 2])
# #   
# #   # THE FOLLWING PART IS ESSENTIAL TO MAKE THE MATRIX CONTAIN THE EXACT SAME IDS AS THE PHENOTYPE FILE
# #   rownames(X.grm) <- nam
# #   colnames(X.grm) <- nam
# #   
# #   attr(X.grm, "rowNames") <- as.factor(rownames(X.grm))
# #   
# #   d <- intersect(ids.object$V2, phenoframe$ANIMALx)  #VECTOR OF IDs shared between phenotype and genotype file
# #   b <- match(phenoframe$ANIMALx, d, nomatch = 20000)
# #   phenoframe <- as.data.frame(subset(phenoframe, b != 20000))  #PHENOTYPE FILE THAT MACHES THE GENOTYPE FILE
# #   
# #   
# #   cn <- match(ids.object[, 2], d, nomatch = 20000)
# #   nam2 <- as.vector(subset(nam, cn != 20000))
# #   
# #   rc <- match(rownames(X.grm), d, nomatch = 20000)
# #   X.grm <- X.grm[rc != 20000, ]  #DELETE REDUNDANT ROWS
# #   cc <- match(colnames(X.grm), d, nomatch = 20000)
# #   X.grm <- X.grm[, rc != 20000]  #DELETE REDUNDANT COLUMNS
# #   X2.grm <- ginv(X.grm)  #CREATE GENERALIZED INVERSE MATRIX
# #   rownames(X2.grm) <- nam2
# #   colnames(X2.grm) <- nam2
# #   attr(X2.grm, "rowNames") <- as.factor(nam2)
# #   
# #   
# #   phenoframe$ANIMALx <- as.factor(phenoframe$ANIMALx)
# #   
# #   x <- list()
# #   x[["grm"]] <- X2.grm
# #   x[["phenoframe"]] <- phenoframe
# #   
# #   x
# # }
# 
# 
# # Provided by Camillo Berenos
# 
# grm.object <- grm.auto
# ids.object <- ids.auto
# head(ids.object)
# 
# id.vector <- ids.object$V2[sample(1:nrow(ids.object), size = 10, replace = F)]
# grm.object <- grm.auto
# ids.object <- ids.auto
# id.vector <- recsumm.f

makeGRM <- function(grm.object, ids.object, id.vector = NULL) {
  
  require(MASS)
  
  elements <- nrow(ids.object)
  X <- diag(elements)
  
  X[upper.tri(X, diag = TRUE)] <- grm.object[,ncol(grm.object)]
  X.grm <- X + t(X) - diag(diag(X))
  nam <- as.factor(ids.object[, 2])
  
  # THE FOLLWING PART IS ESSENTIAL TO MAKE THE MATRIX CONTAIN THE EXACT SAME IDS AS THE PHENOTYPE FILE
  rownames(X.grm) <- nam
  colnames(X.grm) <- nam
  
  if(!is.null(id.vector)) X.grm <- X.grm[which(rownames(X.grm) %in% id.vector),which(colnames(X.grm) %in% id.vector)]
  
  X2.grm <- ginv(X.grm)  #CREATE GENERALIZED INVERSE MATRIX
  rownames(X2.grm) <- rownames(X.grm)
  colnames(X2.grm) <- colnames(X.grm)
  
  attr(X2.grm, "rowNames") <- rownames(X.grm) 

  X2.grm

}
