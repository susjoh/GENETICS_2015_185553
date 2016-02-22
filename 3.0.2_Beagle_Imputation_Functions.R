# 
# gwaa.object <- subHD
# outfile = "test"

GenABEL2beagle <- function(gwaa.object, outfile = NULL){
  
  haplo.obj <- as.character.gwaa.data(gwaa.object)
  #head(haplo.obj)[,1:10]
  
  #~~ Create map
  
  beaglemap <- data.frame(SNP.Name = snp.names(gwaa.object),
                          Position = map(gwaa.object),
                          Ref = refallele(gwaa.object),
                          Eff = effallele(gwaa.object))
  beaglemap$Order <- 1:nrow(beaglemap)
  
  beaglemap <- arrange(beaglemap, Position)
  
  #~~ Create a transposed genotype data frame
  
  temptab <- t(haplo.obj)
  temptab <- cbind(I = rep("M", nrow(temptab)),id = row.names(temptab), temptab)
  
  temptab <- temptab[beaglemap$Order,]
  
  beaglemap <- subset(beaglemap, select = -Order)
  beaglemap$Ref <- gsub(1, "A", beaglemap$Ref)
  beaglemap$Ref <- gsub(2, "C", beaglemap$Ref)
  beaglemap$Eff <- gsub(1, "A", beaglemap$Eff)
  beaglemap$Eff <- gsub(2, "C", beaglemap$Eff)
  
  
  #~~ Convert to Beagle format
  beagletab <- t(apply(temptab, 1, function(x) gsub("/", " ", x, fixed = TRUE)))
  head(beagletab)[,1:5]
  beagletab[is.na(beagletab)] <- "0 0"
  
  colnames(beagletab) <- paste(colnames(beagletab), colnames(beagletab))
  colnames(beagletab)[1:2] <- c("I", "id")
  
  for(i in 3:ncol(beagletab)){
    beagletab[,i] <- gsub("1", "A", beagletab[,i])
    beagletab[,i] <- gsub("2", "C", beagletab[,i])
  }
    
  if(is.null(outfile)){
    list(beagletab = beagletab, beaglemap = beaglemap) 
  } else {
    write.table(beagletab, paste0(outfile, ".txt"), row.names = F, col.names = T, quote = F)
    write.table(beaglemap, paste0(outfile, ".map"), col.names = F, row.names = F, quote = F)
    
  }
  
}



ExtractHaplotypes <- function(input.prefix, results.file, make.ref.alleles = T){
  
  message("Reading and formatting Beagle results file...")
  
  #~~ read in the results file
  
  beagle.res <- read.table(results.file, skip=10, stringsAsFactors = F)
  
  beagle.refalleles <- beagle.res[,4:5]
  
  #~~ transpose data and retain only the haplotype information
  haplotypes.raw <- data.frame(t(beagle.res[,10:ncol(beagle.res)]))
  
  #~~ add ID information and tidy up col names
  
  beagletab <- read.table(paste0(input.prefix, ".txt"), header = T, stringsAsFactors = F)
  beaglemap <- read.table(paste0(input.prefix, ".map"), stringsAsFactors = F)
  
  haplotypes.raw <- cbind(unique(as.vector(dimnames(beagletab)[[2]][seq(3, ncol(beagletab), 2)])),
                          haplotypes.raw)
  names(haplotypes.raw) <- c("ID", as.character(beaglemap[,1]))
  
  #~~ get rid of duplicate ID in each ID string
  
  haplotypes.raw$ID <- as.character(haplotypes.raw$ID)
  haplotypes.raw$ID <- sapply(haplotypes.raw$ID, function(x) strsplit(as.character(x), split = " ")[[1]][1])
  haplotypes.raw$ID <- gsub("X", "", haplotypes.raw$ID)
  haplotypes.raw$ID <- gsub("^\\.", "-", haplotypes.raw$ID)
  
  #~~ convert to character
  
  for(i in 1:ncol(haplotypes.raw)) haplotypes.raw[,i] <- as.character(haplotypes.raw[,i])
  
  #~~ The first two characters (i.e. 0|0) is the phased genotype. 
  #   Retain these and substitute in the reference alleles (if specified).
  
  message("Extracting phased genotypes...")
  
  haplotypes.extract <- data.frame(matrix(NA,
                                          nrow = nrow(haplotypes.raw),
                                          ncol = ncol(haplotypes.raw)))
  
  names(haplotypes.extract) <- names(haplotypes.raw)
  haplotypes.extract$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    haplotypes.extract[,i] <- substr(haplotypes.raw[,i], 1, 3)
    
    if(make.ref.alleles == TRUE){
      haplotypes.extract[,i] <- gsub(0, as.character(beagle.refalleles[i-1,1]), haplotypes.extract[,i])
      haplotypes.extract[,i] <- gsub(1, as.character(beagle.refalleles[i-1,2]), haplotypes.extract[,i])
    }
  }
  
    
  #~~ extract the genotype probabilities
  
  message("Extracting genotype probabilities...")
  
    
  genotype.probs <- data.frame(matrix(NA,
                                      nrow = nrow(haplotypes.raw),
                                      ncol = ncol(haplotypes.raw)))
  names(genotype.probs) <- names(haplotypes.raw)
  genotype.probs$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    genotype.probs[,i] <- sapply(haplotypes.raw[,i], function(x) strsplit(x, split = ":")[[1]][2])
  }
  
  #~~ create a results table for the haplotypes
  
  message("Identifying haplotypes...")
  
  haplo.results <- data.frame(ID = haplotypes.extract$ID,
                              Haplo1 = NA,
                              Haplo2 = NA)
  
  for(i in 1:nrow(haplotypes.extract)){
    haplo.results$Haplo1[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 1, 1), collapse = ""))
    haplo.results$Haplo2[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 3, 3), collapse = ""))
  }
  
  #~~ melt into one column
  require(reshape2)
  haplo.results.1col <- melt(haplo.results, id.vars="ID")
  names(haplo.results.1col) <- c("ID", "HaploID", "Haplo")
  
  for(i in 1:ncol(haplo.results.1col)) haplo.results.1col[,i] <- as.character(haplo.results.1col[,i])
  
  return(list(haplotypes = haplo.results.1col,
              phased.genotypes = haplotypes.extract,
              genotype.probabilities = genotype.probs))
  message("...done.")
}



vector2Dendrogram <- function(haplotype.vector, verbose = TRUE){
  
  #~~ define string similarity function
  
  SeqSimilarity <- function(x1, x2){
    x1 <- as.character(x1)
    x2 <- as.character(x2)
    
    y <- data.frame(x1 = unlist(strsplit(x1[1], split = "")),
                    x2 = unlist(strsplit(x2[1], split = "")),
                    stringsAsFactors = F)
    y$match <- ifelse(y$x1 == y$x2, 1, 0)
    sum(y$match)
  }
  
  haplotype.counts <- data.frame(table(haplotype.vector))
  names(haplotype.counts) <- c("Haplo", "Freq")
  
  #~~ similarity matrix
  
  seq.comp <- data.frame(Seq1temp = rep(as.character(haplotype.counts$Haplo), times = nrow(haplotype.counts)),
                         Seq2temp = rep(as.character(haplotype.counts$Haplo), each  = nrow(haplotype.counts)))
  seq.comp$Seq1 <- apply(seq.comp, MARGIN=1, function (x) (sort(x)[1]))
  seq.comp$Seq2 <- apply(seq.comp[,1:2], MARGIN=1, function (x) (sort(x)[2]))
  
  seq.comp <- unique(subset(seq.comp, select = c(Seq1, Seq2)))
  
  
  seq.comp$similarity.score <- NA
  
  for(i in 1:nrow(seq.comp)){
    seq.comp$similarity.score[i] <- SeqSimilarity(seq.comp$Seq1[i], seq.comp$Seq2[i])
  }
  
  #~~ Create a dendrogram
  
  simi.matrix <- matrix(nrow = nrow(haplotype.counts), ncol=nrow(haplotype.counts))
  simi.matrix[lower.tri(simi.matrix, diag = T)] <- seq.comp$similarity.score
  
  rownames(simi.matrix) <- haplotype.counts$Haplo
  colnames(simi.matrix) <- haplotype.counts$Haplo
  
  simi.matrix
  
  dimnames(simi.matrix)[[1]]
  
  hc = hclust(dist(simi.matrix))
  
  true.vec2 <- which(hc$labels == haplotype.counts$Haplo)
  if(length(true.vec2) == nrow(haplotype.counts)) hc$labels <- paste(haplotype.counts$Freq, haplotype.counts$Haplo, sep = "   ")
  if(length(true.vec2) == nrow(haplotype.counts) & nchar(as.character(haplotype.counts$Haplo[1])) > 30) hc$labels <- haplotype.counts$Freq
  
  library(ggdendro)
  print(ggdendrogram(hc))
  str(hc)
  
  haplotype.counts$Order <- 1:nrow(haplotype.counts)
  
  haplotype.counts$Haplo <- as.character(haplotype.counts$Haplo)
  
  if(verbose == TRUE) list(haplotype.counts = haplotype.counts, hc.object = hc)
}

