# scoring function for differential expression

scoring <- function(data,labels,method="SAM",pcompute="tdist",nperms=1000,memory.limit=TRUE,verbose=TRUE){
  # data: matrix with rows=genes and columns=samples
  # labels: vector or factor of two group labels
  # nperms: number of permutations
  # method: scoring method, currently only "SAM" or "t.test"
  # pcompute : method to compute p-values, either "empirical" from 'nperms'
  #             permutations or "tdist" from t-distribution, or "none" for
  #             no p-value computatation
  # memory limit: is memory limited (less than 1GB)-> slower memory-saving computation
  
  ### 1. Check arguments: ###
  labels<-as.factor(labels)
  if (nlevels(labels)!=2)
    stop("Vector of Labels does not indicate two groups!")
  if (ncol(data)!=length(labels))
    stop("Dimensions of Data and length of Labels-Vector do not match!")
  if (!(pcompute %in% c("empirical","tdist","none")))
    stop("Methode of p-value computation must be \"empirical\",\"tdist\", or \"none\"!\n)")
  binlabels<-as.numeric(labels==levels(labels)[2])  #binary label vector
  ngenes <- nrow(data)
    
  ### 2. Compute observed scores: ###
  tscore <- function(data, group, test="SAM"){
    # group: binary vector of class labels; 1: tumor(tum); 0: control(con)
    # test: one of "SAM" or "t.test"
     revgroup <- 1-group
     ntum <- sum(group)
     ncon <- sum(revgroup)
     meantum <- as.vector(data%*%group/ntum)
     meancon <- as.vector(data%*%revgroup/ncon)
     # Residual Sum of Squares:
     rsstum  <- as.vector((data - matrix(meantum,nrow(data),ncol(data)))^2
                           %*% group)
     rsscon  <- as.vector((data-matrix(meancon,nrow(data),ncol(data)))^2
                           %*% revgroup)
     r     <- meantum - meancon 
     if (test=="SAM"){
       s <- sqrt((1/ncon+1/ntum)/(ncon+ntum-2)*(rsscon+rsstum))
       s0    <- median(s)
       stat  <- r/(s+s0)
     } else if (test=="t.test"){
       # t-Test with equal variances (=SAM without fudge factor)
       s <- sqrt((1/ncon+1/ntum)/(ncon+ntum-2)*(rsscon+rsstum))
       stat  <- r/s
     } # else if (test=="t.test")
     return(stat)
  } #tscore

  # new form to handle permutation matrices (faster than apply, but produces
  #  BIG matrices in between -> often clean workspace, use only 1000 permuations at once
  tscoremat <- function(data,groupmat,test="SAM",memory.limit=TRUE){
    # data: gene expression data with rows=genes and columns=samples
    # groupmat: binary matrix with one column= one vector of classlabels with
    #  1: tumor(tum); 0: control(con)
    nsamples <- ncol(data)
    ngenes <- nrow(data)
    nperms <- ncol(groupmat)    
    onevector <- numeric(nsamples)+1 # for building column-sums
    ntum <- t(groupmat) %*% onevector  # Column-Sums=No.of Tumors
    # normalize group-matrix to column sums=1; divide every 1 by the total number
    #  of ones in that column:
    tummat.norm <- groupmat/matrix(ntum,nrow=nsamples,ncol=nperms,byrow=TRUE)
    # compute means for each gene in tumor and control: 
    meantum <- data %*% tummat.norm
    # compute variance= E(X^2) - (E(X))^2
    vartum <- (data^2) %*% tummat.norm - (meantum^2)
    if (memory.limit) rm(tummat.norm) # clean up
    # the same for control matrix:
    conmat <- 1-groupmat
    ncon <- t(conmat) %*% onevector    # Column-Sums=No.of Controls
    conmat.norm <- conmat/matrix(ncon,nrow=nsamples,ncol=nperms,byrow=TRUE)
    if (memory.limit) rm(onevector,conmat) # clean up
    meancon <- data %*% conmat.norm
    varcon <- (data^2) %*% conmat.norm - (meancon^2)
    r <- meantum - meancon # enumerator of test statistic
    if (memory.limit){rm(data,meantum,meancon,conmat.norm);tempgc <- gc(verbose=FALSE)}
    # compute residual sum of squares: rss=n*var
    rsstum <- matrix(ntum,nrow=ngenes,ncol=nperms,byrow=TRUE)*vartum
    rsscon <- matrix(ncon,nrow=ngenes,ncol=nperms,byrow=TRUE)*varcon
    if (memory.limit){rm(vartum,varcon);tempgc <- gc(verbose=FALSE)}
    s <- sqrt(matrix(((1/ncon+1/ntum)/(ncon+ntum-2)),nrow=ngenes,ncol=nperms,
                     byrow=TRUE)*(rsstum+rsscon))
    if (test=="SAM") {
      s0 <- apply(s,2,median)
      stat <- r/(s+matrix(s0,nrow=ngenes,ncol=nperms,byrow=TRUE))
    } else {
      stat <- r/s
    }# else
    return(stat)
  }# tscoremat
    
  if (verbose) cat("Compute observed test statistics...\n")
  observed.scores <- tscore(data,binlabels,test=method)

  
  ### 3. Compute p-values for observed scores: ###
  if (pcompute=="empirical"){    
    if (verbose) cat("Building permutation matrix...\n")
    perms <- matrix(nrow=nperms,ncol=ncol(data))
    for (i in 1:nperms){  # build matrix of permuted class-labels
      perms[i,] <- sample(binlabels)
    }#for
    if (verbose) cat(paste("Compute",nperms,"permutation test statistics...\n"))
    # compute scores for permuted class-labels:
    if (memory.limit){# use tscoremat but only 250 permutations at once because of memory
      perm.scores <- matrix(nrow=nrow(data),ncol=nperms)
      startcol <- 1
      while (startcol < nperms){
        endcol <- min(nperms,startcol+249)      
        perm.scores[,startcol:endcol] <- tscoremat(data,t(perms)[,startcol:endcol],test=method)
        if (verbose) cat(endcol,"...")
        startcol <- startcol+250
      } # while (startcol < nperms)
    } else 
       perm.scores <- tscoremat(data,t(perms),test=method)
    
    if (verbose) cat("\nCompute empirical p-values...\n")
    # how many permutation scores are absolutely greater or equal
    #  than the respective observed score?
    pvalues <- rowSums(abs(perm.scores) >= matrix(abs(observed.scores),nrow=ngenes,ncol=nperms,byrow=FALSE))/nperms

  } else if (pcompute=="tdist"){
    # Auxiliary Function to compute p-values for given t-distribution;
    # can adjust for multiple testing with method of Holm et al.
    pval <- function(tscores,t.df,adjust=TRUE){
      rawp <- sapply(tscores,pt,df=t.df,lower.tail=FALSE)
      if (adjust){
        rankp <- rank(rawp,ties.method="first")
        sortedp <- sort(rawp,method="quick")
        multp <- sortedp * rev(seq(length(tscores)))
        multp[multp>1] <- 1
        adjustedp <- multp[rankp]
      } else {
        adjustedp <- rawp
      } #else
      return(adjustedp)
    } #pval    
    if (verbose) cat("Computing p-values from t-distribution.\n")
    pvalues <- pval(observed.scores,t.df=ncol(data)-2)
  } else if (pcompute=="none") {
    pvalues <- rep(NA,length(observed.scores))
  } # ifelse pcompute

  ### 4. prepare and return result: ###
  if (verbose) cat("Compute quantiles of empirical distributions...")
  # return quantiles of permuted score matrix as borders for expected scores:
  if (pcompute=="empirical") {
    expected.borders <- apply(perm.scores,1,quantile,probs=c(0.025,0.975))
    expected.lower <- expected.borders[1,]
    expected.upper <- expected.borders[2,]
    if (!(is.null(rownames(data))))
      names(expected.lower) =  names(expected.upper) <- rownames(data)
    else 
      names(expected.lower) =  names(expected.upper) <- seq(nrow(data))
  } else {
    expected.lower = expected.upper <- NULL
  } # ifelse (pcompute=="empirical")
  
  if (!(is.null(rownames(data))))
    names(observed.scores)=names(pvalues) <- rownames(data)
  else 
    names(observed.scores)=names(pvalues) <- seq(nrow(data))  

  result <- list(observed=observed.scores, pvalues=pvalues,
                 expected.lower=expected.lower,expected.upper=expected.upper)
  if (verbose) cat("Done.\n")
  return(result)  
} #scoring
