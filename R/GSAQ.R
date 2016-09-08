###################################Dependent Packages#####################################
library(stats)
##########################################################################################
#      Gene set enrichment analysis with Quantitative Trait Loci data                                                                                  #
#                                                                                        #
##########################################################################################

#############################Selection of informative genes###############################

GeneSelect <- function(x, y, s, method=c("t-score", "F-score", "MRMR", "BootMRMR"))

{
  this.call = match.call()
  if ((!class(x)=="data.frame")) {
    warning("x must be a data frame and rows as gene names")
    }   
  if ((!class(y)=="numeric")) {
      warning("y must be a vector of 1/-1's for two class problems")
      }   
   if (!length(y)==ncol(x))
   {
     warning("Number of samples in x must have same number of sample labels in y")
   }
  if (!(s > 0 & s < nrow(x)))
  {
    warning ("size of the selected gene set must be numeric and less than total number of genes")
  }
  if ((!class(method)=="character"))
  {
    warning ("method must be the name of the method to be used for gene selection")
  }
  
  cls <- as.vector(y) ###class information###
  genenames <- rownames(x)
  g <- as.matrix(x)
  n <- nrow(g)               ###number of genes####
  M <- ncol(x)           ###number of samples###
  method <- match.arg(method)
  if (method=="t-score")
  {
    cls <- as.numeric(y)         
    ranking <- vector(length=n)
    idx <- which(cls==1)  # indexing of positive samples
    idy <- which(cls==-1) # indexing of negative samples
    #print(idx)
    B=vector(mode="numeric", n)
    for(i in 1:n){
      t.mes <-t.test(g[i, idx], g[i, idy], alternative="two.sided")$statistic  #####F-Score    
      B[i] <- t.mes
    }
    ranking = sort(-B, index.return = TRUE)$ix
    temp <- ranking[1:s]
    selectgenes <- genenames[temp]
    out1 <- list(selectgenes=selectgenes, method=method)
    class(out1) <- c("Informative geneset")
    #selectgenes
    out1
  }
  if(method=="F-score")
  {
    cls <- as.numeric(y)  
    idx <- which(cls==1)  # indexing of positive samples
    idy <- which(cls==-1) # indexing of negative samples
    B=vector(mode="numeric", n)
    for(i in 1:n){
      f.mes <-(((mean(g[i, idx])-mean(g[i, ]))^2)+ ((mean(g[i, idy])-mean(g[i, ]))^2))/(var(g[i, idx])+var(g[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    ranking <- sort(-B, index.return = TRUE)$ix 
    temp <- ranking[1:s]
    selectgenes <- genenames[temp]
    out1 <- list(selectgenes=selectgenes, method=method)
    class(out1) <- c("Informative geneset")
    #selectgenes
    out1
  }
  
  if (method=="MRMR")
  {
    cls <- as.numeric(y) ###class information###
    GeneRankedList <- vector(length=n)
    qsi <- as.vector((apply(abs(cor(t(g), method="pearson",use="p")-diag(n)), 1, sum))/(n-1))
    idx <- which(cls==1)  # indexing of positive samples
    idy <- which(cls==-1) # indexing of negative samples
    B <- vector(mode="numeric", n)
    for(i in 1:n){
      f.mes <-(((mean(g[i, idx])-mean(g[i, ]))^2)+ ((mean(g[i, idy])-mean(g[i, ]))^2))/(var(g[i, idx])+var(g[i, idy]))  #####F-Score    
      B[i] <- f.mes
    }
    rsi <- abs(B)
    rankingCriteria <- rsi/qsi
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    temp <- GeneRankedList[1:s]
    selectgenes <- genenames[temp]
    out1 <- list(selectgenes=selectgenes, method=method)
    class(out1) <- c("Informative geneset")
    #selectgenes
    out1
  }
  
  if (method=="BootMRMR")
  {
    cls <- as.numeric(y)               ### class information###
    genes <- rownames(x)
    g <- as.matrix(x)
    n <- nrow(g)                       ### number of genes####
    M <- ncol(x)                       ### number of samples###
    S <- 85                            ### number of bootstraps ###
    Q <- 0.5                           ### Quartile value under null hypothesis###
    GeneRankedList <- vector(length=n)
    M1 <- matrix(0, n, S) 
    ##if(missing(s)) 
    for (j in 1:S) {
      samp <- sample(M, M, replace=TRUE) ###select bootstrap sample #####
      x1 <- g[, samp]
      y1 <- cls[samp] 
      qsi <- as.vector((apply(abs(cor(t(x1), method="pearson",use="p")-diag(n)), 1, sum))/(n-1))
      idx <- which(y1==1)                # indexing of positive samples
      idy <- which(y1==-1)               # indexing of negative samples
      B=vector(mode="numeric", n)
      for(i in 1:nrow(x1)){
        f.mes <-(((mean(x1[i, idx])-mean(x1[i, ]))^2)+ ((mean(x1[i, idy])-mean(x1[i, ]))^2))/(var(x1[i, idx])+var(x1[i, idy]))  #####F-Score    
        B[i] <- f.mes
      }
      rsi <- abs(B)
      rankingCriteria <- rsi/qsi
      GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
      rankvalue <- sort(GeneRankedList, index.return=TRUE)$ix
      rankscore <- (n+1-rankvalue)/(n)
      M1[,j] <- as.vector(rankscore)
    }
    rankscore <- as.matrix(M1)          # Raw scores obtained from SVM-RFE
    mu <- Q                                    # value under null hypothesis
    
    R <- rankscore - mu                       # Transformed score of MRMR under H0
    sam <- nrow (R)                           # number of genes
    pval.vec <- vector(mode="numeric", length=nrow(rankscore))
    for (i in 1:sam) {
      z <- R[i,]
      z <- z[z != 0]
      n11 <- length(z)
      r <- rank(abs(z))
      tplus <- sum(r[z > 0])
      etplus <- n11 * (n11 + 1) / 4
      vtplus <- n11 * (n11 + 1) * (2 * n11 + 1) / 24
      p.value=pnorm(tplus, etplus, sqrt(vtplus), lower.tail=FALSE)
      pval.vec[i]=p.value
    }
    w11 <- as.vector(pval.vec)
    gene.id <- sort(w11, index.return=TRUE)$ix
    temp <- gene.id [1:s]
    selectgenes <- genes[temp]
    out1 <- list(selectgenes=selectgenes, method=method)
    class(out1) <- "Informative geneset"
    out1
  }
  if(missing(method))
  {
    method=="t-score"
  }
  return(out1)
  }
################## Computation of QTL hits in the selected gene set #####################

qtlhit <- function(geneset, genelist, qtl)
{
  this.call = match.call()
  
  if ((!class(geneset)=="character")) {
    warning("gene set must be a vector of gene names")
  }
  if((!class(qtl)=="data.frame")) {
    warning("QTL must be a dataframe with row names as QTL ids and three columns as chromosome number, start and stop positions")
  }
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if(length(geneset) > nrow(genelist)){
    warning("geneset must be a sub set of genelist")
  }
  geneset <- as.vector(geneset)
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
  } 
  tempp <- temp[!is.na(temp)] 
  select.gen <- genelist[tempp,] 
  chr.gen <- select.gen[,1]
  chr.qtl <- qtl[,1]
  gene.start <- as.vector(select.gen[,2])
  gene.stop <- as.vector(select.gen[,3])
  qtl.start <- as.vector(qtl[,2])
  qtl.stop <- as.vector(qtl[,3])
  m <- nrow(qtl) #######number of QTLs####
  n <- nrow (geneset) #######number of genes####
  selgenenam <- rownames(select.gen)
  
  Mat <- NULL
  
  for (i in 1:length(selgenenam))
  {
    id <- as.vector(which(chr.qtl==chr.gen[i]))
    if(length(id)>0)
    {
      a <- matrix("NULL",1, length(id))
      for (j in 1:length(id))
      {
        a=ifelse(gene.start[i] >= qtl.start[id[j]] & gene.stop[i]<= qtl.stop[id[j]], 1, 0)
        Mat <- rbind(Mat,a)      
      }     
    } 
  }  
  Mat <- as.vector(Mat)
  qtlhits <- sum(Mat==1)
  if (sum(Mat==1)==0)
  {
    warning(" Selected geneset has no qtlhits. Consider another geneset with larger size")
  }
  class( qtlhits) <- "Total QTLhits in Geneset"
  return(qtlhits)
}

########################Gene Set Validation with QTL data without gene sampling###########

GSVQ <- function (geneset, genelist, qtl) 
{
  this.call = match.call()
  
  if ((!class(geneset)=="character")) {
    warning("gene set must be a vector of gene names")
  }
  if((!class(qtl)=="data.frame")) {
    warning("QTL must be a dataframe with row names as QTL ids and three columns as chromosome number, start and stop positions")
  }
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if(length(geneset) > nrow(genelist)){
    warning("geneset must be a sub set of genelist")
  }
  
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
  }
  
  totqtlhit <- function(genelist, qtl){
    genename1 <- rownames(genelist)
    tothits <- qtlhit(genename1, genelist, qtl)
    return(tothits)
  }
  m1 <- totqtlhit (genelist, qtl)
  no.qtl <- qtlhit(geneset, genelist, qtl)
  pval1 <- phyper(no.qtl, m1, length(genenames)-m1, length(geneset), lower.tail=F)
  res <- list(no.qtl=no.qtl, pval1=pval1 )
  class(res) <- c("Total QTL Hits", "P-value")
  res
} 

#################Total QTL hits found in the whole gene list#############################

totqtlhit <- function (genelist, qtl)
{
  if((!class(qtl)=="data.frame")) {
    warning("QTL must be a dataframe with row names as QTL ids and three columns as chromosome number, start and stop positions")
  }
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  genename1 <- rownames(genelist)
  tothits <- qtlhit(genename1, genelist, qtl)
  class(tothits)="Total QTLhits in genelist"
  return(tothits)
}
  
######################################Gene Set Analysis with QTL data ###################

GSAQ <- function (geneset, genelist, qtl, SampleSize, K, method=c("meanp", "sump", "logit", "sumz", "logp")) 
{
  this.call = match.call()
  
  if ((!class(geneset)=="character")) {
    warning("gene set must be a vector of gene names")
  }
  if((!class(qtl)=="data.frame")) {
    warning("QTL must be a dataframe with row names as QTL ids and three columns as chromosome number, start and stop positions")
  }
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if(length(geneset) > nrow(genelist)){
    warning("geneset must be a sub set of genelist")
  }
  if(length(geneset) < SampleSize){
    stop("sample size must be less than the size of the selected gene set")
  }
  method <- match.arg(method)
  if(K <=0){
    stop("Number of samples drawn must be greater than zero")
  }
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the genelist. Make sure\that gene identifiers in dataset are Gene symbols")
  }
  
  totqtlhit <- function (genelist, qtl)
    {
    genename1 <- rownames(genelist)
    tothits <- qtlhit(genename1, genelist, qtl)
    return(tothits)
    }
  m1 <- totqtlhit (genelist, qtl)
  qtl1 <- function(geneset1, genelist, qtl)
    {
    geneset1 <- as.vector(geneset1)
    qtl.start <- as.vector(qtl[,2])
    qtl.stop <- as.vector(qtl[,3])
    m <- nrow(qtl) #######number of QTLs####
    n <- nrow (geneset1) #######number of genes####
    genenames <- rownames(genelist) 
    temp1 <- match(unlist(geneset1), genenames)
    tempp <- temp1[!is.na(temp1)] 
    select.gen <- genelist[tempp,] 
    chr.gen <- select.gen[,1]
    chr.qtl <- qtl[,1]
    gene.start <- as.vector(select.gen[,2])
    gene.stop <- as.vector(select.gen[,3])
    chr.slt <- as.numeric(names(table(chr.qtl)))
    selgenenam <- rownames(select.gen)
    
    Mat <- NULL
    
    for (i in 1:length(selgenenam))
    {
      id <- as.vector(which(chr.qtl==chr.gen[i]))
      if(length(id)>0)
      {
        a <- matrix("NULL",1, length(id))
        for (j in 1:length(id))
        {
          if(gene.start[i] >= qtl.start[id[j]] & gene.stop[i]<= qtl.stop[id[j]])
          {
            
            a <- rownames(select.gen[i,])
            Mat <- rbind(Mat,a)
          }       
        }     
      } 
    }  
    Mat <- as.vector(Mat)
    #print(Mat)
    no.qtl <- length(Mat)
    #print(no.qtl)
    return(no.qtl)
  }
  size <- vector(length=K)
  for( i in 1:length(size)){
    samp <- sample(1:length(geneset), SampleSize, replace=F)
    geneset1 <- as.vector(geneset[samp])
    hit1 <- qtl1(geneset1, genelist, qtl)
    pval <- phyper(hit1, m1, length(genenames)-m1, length(samp), lower.tail=F)
    size[i]<- pval
  }
  
  if(method=="logit")
      {
      p <- size
      k1 <- (p >= 0) & (p < 1)
      psum <- sum(log(p[k1]/(1 - p[k1])))
      k <- length(p[k1])
      if (sum(1L * k1) < 2) 
        stop("Must have at least two valid p values")
      mult <- -1/sqrt(k * pi^2 * (5 * k + 2)/(3 * (5 * k + 4)))
      
      t <- mult * psum
      df <- (5 * k + 4)
      res <- list(t = t, df = df, p = pt(t, df, lower.tail = FALSE), method=method)
      class(res) <- c("logitp","df", "finalp")
      res
    }
  
  ###### Mean Method####
  if(method=="meanp")
  {
    p <- size
    k1 <- (p >= 0) & (p <= 1)
    pi <- mean(p[k1])
    k <- length(p[k1])
    if (sum(1L * k1) < 4) 
      stop("Must have at least four valid p values")
    z <- (0.5 - pi) * sqrt(12 * k)
    res <- list(z = z, p = pnorm(z, lower.tail = FALSE), method=method)
    class(res) <- c("meanp", "finalp")
    res
  }
  
 if (method=="sumz")
 {
    p <- size
     k1 <- (p >= 0) & (p <=1)
     if (sum(1L * k1) < 2) 
       stop("Must have at least two valid p values")
        weights <- rep(1, length(p))
     if (sum(1L * k1) != length(p)) {
       warning("Some samples are omitted")
       omitw <- weights[!k1]
       if (sum(1L * omitw) > 0) 
         warning("Weights omitted too")
     }
     zp <- (qnorm(p[k1], lower.tail = FALSE) %*% weights[k1])/sqrt(sum(weights[k1]^2))
     res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE), method=method) 
                
     class(res) <- c("sumz", "finalp")
     res
   } 
 
 
 if (method=="sump")
   {
   p <- size  
   k1 <- (p >= 0) & (p <= 1)
     sigmap <- sum(p[k1])
     k <- length(p[k1])
     if (sum(1L * k1) < 2) 
       stop("Must have at least two valid p values")
     conservativep <- exp(k * log(sigmap) - lgamma(k + 1))
     nterm <- floor(sigmap) + 1
     denom <- lfactorial(k)
     psum <- 0
     terms <- vector("numeric", nterm)
     for (i in 1:nterm) {
       terms[i] <- lchoose(k, i - 1) + k * log(sigmap - i + 
                                                 1) - denom
       pm <- 2 * (i%%2) - 1
       psum <- psum + pm * exp(terms[i])
     }
     if (k != length(p)) {
       warning("Some samples omitted")
     }
     if (sigmap > 20) 
     {
       warning
     }("Likely to be unreliable, check with another method")
     res <- list(p = psum, conservativep = conservativep, method=method)
     class(res) <- c("sump", "finalp")
     res
   }
 
 if (method=="logp")
 {
   p <- size
   k1 <- (p > 0) & (p <= 1)
   lnp <- log(p[k1])
   chisq <- (-2) * sum(lnp)
   df <- 2 * length(lnp)
   if (sum(1L * k1) < 2) 
     stop("Must have at least two valid p values")
   if (length(lnp) != length(p)) {
     warning("Some samples omitted")
   }
   res <- list(chisq = chisq, df = df, p = pchisq(chisq, df, 
                lower.tail = FALSE), method=method)
   class(res) <- c("sumlog", "df", "finalp")
   res
 }
 if(missing(method))
 {
   method=="logp"
 }
 return(res)
 }

 ########### Chromosomal distribution of selected gene set
 
 genedist <- function (geneset, genelist, plot =TRUE)
 {
   this.call = match.call()
   
   if ((!class(geneset)=="character")) {
     warning("gene set must be a vector of gene names")
   }
   
   if(!class(genelist)=="data.frame") {
     warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
   }
   if(length(geneset) > nrow(genelist)){
     warning("geneset must be a sub set of genelist")
 }
 
 genenames <- rownames(genelist)
 temp <- match(unlist(geneset), genenames)
 if (sum(is.na(temp)) == length(temp)) {
   stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
 }
 if (sum(!is.na(temp))/length(temp) < 0.05) {
   stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
 }
 tempp <- temp[!is.na(temp)]
 selectgene <- genelist[tempp,]
 selectgennam <- genenames[tempp]
 
 if(plot==TRUE){
   g <- selectgene[,1]
   gg <- table(g)
   chr.nam <- as.numeric(names(gg))
   plot1 <- hist(g, breaks=max(chr.nam), freq=F, prob=T, right=FALSE, 
        col=rainbow(length(chr.nam)), main="Chromosomal distribution of selected genes", 
        xlim=c(1, max(chr.nam)), xlab="Chromosome number")
   lines(density(g)) 
 }
 out1 <- as.data.frame(cbind(selectgennam, selectgene))
 rownames (out1) <- NULL
 colnames(out1) <- c("Gene ID", "Chr. No.", "Start Position (bp)", "End Position (bp)")
 #class(out1) <- "Genomic Positions of Selected Genes"
 return(out1)
 }
 
###############################Positions of genes in respective QTLs######################

geneqtl <- function (geneset, genelist, qtl)
{
  this.call = match.call()
  
  if ((!class(geneset)=="character")) {
    warning("gene set must be a vector of gene names")
  }
  
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if(length(geneset) > nrow(genelist)){
    warning("geneset must be a sub set of genelist")
  }
  
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
  }
  
  qtl.start <- as.vector(qtl[,2])
  qtl.stop <- as.vector(qtl[,3])
  m <- nrow(qtl) #######number of QTLs####
  n <- nrow (geneset) #######number of genes####
  geneset <- as.vector(geneset)
  genenames <- rownames(genelist) 
  temp1 <- match(unlist(geneset), genenames)
  tempp <- temp1[!is.na(temp1)] 
  select.gen <- genelist[tempp,] 
  chr.gen <- select.gen[,1]
  chr.qtl <- qtl[,1]
  gene.start <- as.vector(select.gen[,2])
  gene.stop <- as.vector(select.gen[,3])
  chr.slt <- as.numeric(names(table(chr.qtl)))
  selgenenam <- rownames(select.gen)
  
  Mat <- NULL
  Mat1 <- NULL
  
  for (i in 1:length(selgenenam))
  {
   id <- as.vector(which(chr.qtl==chr.gen[i]))
   if(length(id)>0)
   {
     a <- matrix("NULL",1, length(id))
     b <- matrix("NULL", 1, length(id))
     for (j in 1:length(id))
     {
      if(gene.start[i] >= qtl.start[id[j]] & gene.stop[i]<= qtl.stop[id[j]])
      {
     #print(select.gen[i,])
     #print(rownames(qtl[id[j],]))
        a <- cbind(rownames(select.gen[i,]), rownames(qtl[id[j],]))
        b <- cbind(rownames(select.gen[i,]), select.gen[i,], rownames(qtl[id[j],]), qtl[id[j],])
       #print(a) 
       #print(b)
       Mat <- rbind(Mat,a)
       Mat1 <- rbind(Mat1, b)
      }       
     }     
   } 
  }  
  colnames (Mat) <- c("Gene ID", "QTL ID")
  Mat <- as.data.frame (Mat)
  rownames(Mat1) <- NULL
  colnames(Mat1) <- c("Gene ID", "Chr. No.", "Start (bp)", "End (bp)", "QTL ID", "Chr. No", "Start (bp)", "End (bp)")
  
  res <- list(Mat=Mat, Mat1=Mat1)
  class(res) <- c("List of qtlhit genes", "Genomic positions of qtlhit genes")
  res
}

##############################QTL wise distribution of selected genes##################
qtldist <- function(geneset, genelist, qtl, plot=TRUE)
{
  this.call = match.call()
  
  if ((!class(geneset)=="character")) {
    warning("gene set must be a vector of gene names")
  }
  
  if(!class(genelist)=="data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if(length(geneset) > nrow(genelist)){
    warning("geneset must be a sub set of genelist")
  }
  
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
  }
  
  qtl.start <- as.vector(qtl[,2])
  qtl.stop <- as.vector(qtl[,3])
  m <- nrow(qtl) #######number of QTLs####
  n <- nrow (geneset) #######number of genes####
  geneset <- as.vector(geneset)
  genenames <- rownames(genelist) 
  temp1 <- match(unlist(geneset), genenames)
  tempp <- temp1[!is.na(temp1)] 
  select.gen <- genelist[tempp,] 
  chr.gen <- select.gen[,1]
  chr.qtl <- qtl[,1]
  gene.start <- as.vector(select.gen[,2])
  gene.stop <- as.vector(select.gen[,3])
  chr.slt <- as.numeric(names(table(chr.qtl)))
  selgenenam <- rownames(select.gen)
  
  Mat <- NULL
  
  for (i in 1:length(selgenenam))
  {
    id <- as.vector(which(chr.qtl==chr.gen[i]))
    if(length(id)>0)
    {
      a <- matrix("NULL",1, length(id))
      for (j in 1:length(id))
      {
        if(gene.start[i] >= qtl.start[id[j]] & gene.stop[i]<= qtl.stop[id[j]])
        {

          a <- cbind(rownames(select.gen[i,]), rownames(qtl[id[j],]))
          Mat <- rbind(Mat,a)
        }       
      }     
    } 
  }  
  colnames (Mat) <- c("Gene ID", "QTL ID")
  Mat <- as.data.frame (Mat)
  
  if (plot==TRUE)
  {
    qtl11 <- table(Mat [,2])
    QTLID <- names (qtl11)
    cc <- as.factor(Mat[,2])
    plot <- barplot(table(cc), main="QTL Wise Distribution QTLhit Genes",
                    xlab = "QTL IDs", ylab = "Number of QTLhit Genes",ylim = c(0, max(as.vector(qtl11))),
                    col=rainbow(length(QTLID)), names=as.vector(QTLID), space=1, las=2)
  }
  class(qtl11) <- "QTL wise distribution of selected genes"
  return(qtl11)
  
  
}
#################################Ends here################################################

