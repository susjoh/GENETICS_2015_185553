
CumuPos<-function(scan.gwaa.object) {
  

  
  #~~ Extract results and calculate cumulative positions
  gwa_res <- results(scan.gwaa.object)
  
  gwa_res$Chromosome <- as.character(gwa_res$Chromosome)
  
  if("X" %in% gwa_res$Chromosome){
    chrvec <- unique(gwa_res$Chromosome)
    
    if(length(chrvec == 1)){
      print("Assuming X chromosome is chromosome 27 - this function will be fixed in future")
      gwa_res$Chromosome[which(gwa_res$Chromosome == "X")] <- 27
    }
      
    if(length(chrvec > 1)){
      chrvec <- as.numeric(chrvec[-which(chrvec == "X")])
      gwa_res$Chromosome[which(gwa_res$Chromosome == "X")] <- max(chrvec) + 1
    }
    
  }
  
  gwa_res$Chromosome <- as.numeric(gwa_res$Chromosome)
  
  gwa_res <- gwa_res[with(gwa_res, order(Chromosome,Position)), ]
  
  gwa_res$Diff <- c(0,diff(gwa_res$Position))
  gwa_res$Diff[gwa_res$Diff < 0] <- 1

  gwa_res$Cumu <- cumsum(gwa_res$Diff)
  
  gwa_res$Cumu2 <- gwa_res$Cumu + (25000000 * gwa_res$Chromosome)
  
  
  #~~ Plot the uncorrected p-values
  return(gwa_res)
}



FullGwasPlot<-function(cumu.object, corrected = FALSE, bonf = F, lambda = NULL) {
  
  require(ggplot2)
  
  #chrinfo <- data.frame(Chromosome = 0, Start = 0, Stop = 0)
  chrinfo <- NULL
  
  if(bonf == F) bonf = 0.05/nrow(cumu.object)
  
  for(i in unique(cumu.object$Chromosome)){
    
    temp1 <- subset(cumu.object, Chromosome == i)
    
    temp2 <- data.frame(Chromosome = i,
                        Start = temp1[1,"Cumu"],
                        Stop = temp1[nrow(temp1),"Cumu"])
    
    chrinfo <- rbind(chrinfo, temp2)
  }
  
  chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
  
  colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$Chromosome)))
  colourscale <- colourscale[1:(length(colourscale)/2)]
  
  
  if(corrected){  
    
    if(!is.null(lambda)){
      cumu.object$Pc1df <- 1-pchisq(cumu.object$chi2.1df/lambda, 1)
      warning("lambda adjusted in plot only")
    }
    
    ggplot(cumu.object, aes(Cumu, -log10(Pc1df), col = factor(Chromosome))) +
      geom_point(size = 3.5,alpha = 0.4) +
      geom_hline(yintercept = -log10(bonf),linetype = 2, alpha = 0.6, size = 1) +
      scale_colour_manual(values = colourscale) +
      theme(legend.position = "none") +
      theme(axis.text.x  = element_text (size = 16, vjust = 0),
            axis.text.y  = element_text (size = 14, hjust = 1.3),
            strip.text.x = element_text (size = 16, vjust = 0.7),
            axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
            axis.title.x = element_text (size = 16, vjust = 0.2),
            strip.background = element_blank()) +
      scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chromosome) +
      labs(x = "Chromosome", y = "-log10 P")
    
  } else {
    
    ggplot(cumu.object, aes(Cumu,-log10(P1df), col = factor(Chromosome))) +
      geom_point(size = 3.5, alpha = 0.4) +
      geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
      scale_colour_manual(values = colourscale) +
      theme(legend.position="none") +
      theme(axis.text.x  = element_text (size = 16, vjust = 0),
            axis.text.y  = element_text (size = 14, hjust = 1.3),
            strip.text.x = element_text (size = 16, vjust = 0.7),
            axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
            axis.title.x = element_text (size = 16, vjust = 0.2),
            strip.background = element_blank()) +
      scale_x_continuous(breaks=chrinfo$Mid,labels=chrinfo$Chromosome) +
      labs(x="Chromosome",y="-log10 P")
  }
  
  
}



FullPpPlot<-function(scan.gwaa.object = NULL, cumu.object, corrected = FALSE, lambda = NULL) {
  
  require(ggplot2)
  
  if(corrected){
    
    if(!is.null(lambda)){
      cumu.object$Pc1df <- 1-pchisq(cumu.object$chi2.1df/lambda, 1)
      warning("lambda adjusted in plot only")
    }
    
    gwa_res.null<-data.frame(obs=sort(-log10(cumu.object$Pc1df)),
                             exp=sort(-log10(seq(1/nrow(cumu.object),1,1/nrow(cumu.object)))))
    
  } else {
    
    gwa_res.null<-data.frame(obs=sort(-log10(cumu.object$P1df)),
                             exp=sort(-log10(seq(1/nrow(cumu.object),1,1/nrow(cumu.object)))))
    
  }
  
  ggplot(gwa_res.null,aes(x=exp,y=obs)) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          strip.background = element_blank()) +
    labs(x="Expected -log10 P",y="Observed -log10 P")
  
} 

mdsPlot <- function(mdsframe, trait){
  mdsframe$trait <- mdsframe[,trait]
  
  require(ggplot2)
  ggplot(mdsframe, aes(X1, X2, col = factor(trait))) + 
    geom_point(size = 3, alpha = 0.7) +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2)) +
    labs(title = paste("Clustered by",trait)) +
    scale_color_brewer(palette = "Set1")
}


FullSummary <- function(scan.gwaa.object, bonf = F) {
  
  cumutemp <- CumuPos(scan.gwaa.object)
  
  if(bonf == F) bonf = 0.05/nrow(cumutemp)
  
  
  print(summary(scan.gwaa.object))
  #print(lambda(scan.gwaa.object))
  
  p1 <- FullGwasPlot (cumutemp) + labs(title = "Uncorrected GWAS")
  p2 <- FullPpPlot (scan.gwaa.object,  cumutemp) + labs(title = paste("lambda =",lambda(scan.gwaa.object)$estimate))
  
  p3 <- FullGwasPlot (cumutemp, corrected = TRUE) + labs(title = "Corrected GWAS")
  p4 <- FullPpPlot (scan.gwaa.object, cumutemp,  corrected = TRUE) + labs(title = "lambda = 1")
  
  multiplot(p1, p2, p3, p4, cols = 2) 
  
}

