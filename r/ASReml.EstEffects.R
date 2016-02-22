
#model <- teno.model

ASReml.EstEffects <- function(model){
  x <- summary(model)$varcomp
  totvar   <- sum(x[,1])
  x$Effect <- x$gamma/totvar
  x$SE <- NA
  
  object <- model
  pframe <- as.list(object$gammas) 
  
  denominator <- "1 "
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    denominator <- paste(denominator, "+ `", effname, "` ", sep = "")
  }
  
  denominator <- paste("(", denominator, ")", sep = "")
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    transform <- eval(parse(text = paste("`placeholder` ~ `", effname, "`/", denominator, sep = "")))
    
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe) 
    X <- as.vector(attr(tvalue, "gradient")) 
    tname <- if (length(transform) == 3) {
      transform[[2]] 
    } else "" 
    V <- object$ai 
    n <- length(pframe) 
    i <- rep(1:n, 1:n) 
    if (length(V) != length(i)) 
      stop("vcov matrix incompatible with\nnumber of variance components") 
    j <- sequence(1:n) 
    k <- 1 + (i > j) 
    se <- sqrt(sum(V * X[i] * X[j] * k))
    x[effname,"SE"] <- se
  }
  
  #~~ NEED TO FIX THIS ADDITION FROM JISCA
  
if(nrow(summary(model)$varcomp) > 2){
  VC <- summary(model)$varcomp
  PrV <-  setNames(VC[,2] / sum(VC[,2]), rownames(VC))
  n <- nrow(VC)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  vcov <- matrix(NA,n,n)
  vcov[cbind(i,j)] <- model$ai
  vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
  vcovX <- vcov[-n,-n]   # exclude residual variance; irrelevant for proportions
  PrV.SE <- numeric()
  for (k in 1:(n-1)) {
    PrV.SE[k] <- sqrt(1/sum(VC$gamma)^4 *
                        (sum(VC$gamma[-k])^2 * diag(vcovX)[k] +
                           VC$gamma[k]^2 * sum(vcovX[-k,-k]) - # diagonal & off-diagonal
                           2*VC$gamma[k] * sum(VC$gamma[-k]) * sum(vcovX[k,-k])))
  }
  PrV.SE[n] <- sqrt(1/sum(VC$gamma)^4 * sum(vcovX))
  
  #~~~
  
  x[nrow(x), ncol(x)] <- PrV.SE[n]
  }
  
  x
}


