
crossvalParam <- function(locations,scores,kernel,kernelparams=NULL,ncross=10,verbose=TRUE){  
  # do cross-validation:
  nscores <- length(scores)
  setidx <- split(sample(1:nscores),rep(1:ncross,length.out=nscores))
  totalresid <- 0
  for (i in 1:length(setidx)){
    if (verbose) cat("Iteration",i,"...")
    trainScores <- scores[-setidx[[i]]]
    testScores <- scores[setidx[[i]]]
    trainlocs <- locations[-setidx[[i]]]
    testlocs <- locations[setidx[[i]]]
    kernelweights <- kernelmatrix(testlocs,trainlocs,kernel,kernelparams)
    slideScores <- kernelize(trainScores,kernelweights)
    if (verbose) cat("Computing residuals for test scores...\n")
    resid <- squaredDist(testScores,slideScores)
    totalresid <- totalresid + resid
  } # for (i in 1:length(setidx))
  avgresid <- totalresid/length(setidx)
  #avgavgsq <- totalavgsq/nscores # normalize by number of genes on chromosome
  return(avgresid)
} # crossvalParam

evalParams <- function(locations,scores,kernel,kernelparams,paramMultipliers=2^(-4:4),ncross=10,verbose=TRUE){
  parameterName <- names(kernelparams)[1]
  tryparams <- c()
  tryresid <- c()
  cparams <- kernelparams
  # evaluate multipliers of parameter for cross-validation performance:
  for (thisMultiplier in paramMultipliers){
    if (identical(kernel,kNN))
      cparams[[1]] <- ceiling(kernelparams[[1]]*thisMultiplier)
    else
      cparams[[1]] <- kernelparams[[1]]*thisMultiplier
    thisParam <- cparams[[1]][1]
    if (verbose) cat("\nEvaluating parameter",parameterName,"=",thisParam,"...\n")
    multiplierError <- crossvalParam(locations=locations,scores=scores,kernel=kernel,kernelparams=cparams,ncross=ncross,verbose=verbose)
    tryparams <- c(tryparams,thisParam)
    tryresid <- c(tryresid,multiplierError)
  } # for (thisMultiplier in paramMultipliers)
  bestParam <- list(x=tryparams[which.min(tryresid)])
  names(bestParam)[1] <- parameterName
  result <- list(parameterName=tryparams,avgResid=tryresid,multiplier=paramMultipliers,best=bestParam)
  names(result)[1] <- parameterName
  class(result) <- "MACATevP"
  return(result)
} # evalParams

evaluateParameters <- function(data,class,chromosome,kernel,kernelparams=NULL,paramMultipliers=2^(-4:4),subset=NULL,newlabels=NULL,ncross=10,verbose=TRUE){
  ### check arguments: ###
  stopifnot(!is.null(data$expr),!is.null(data$labels),
            !is.null(data$geneName),!is.null(data$geneLocation),
            !is.null(data$chromosome),!is.null(data$chip))
  if(is.null(subset))
    subset <- seq(ncol(data$expr))
  if(!is.null(newlabels)){
    if(length(newlabels)!=length(subset))
      stop("Length of new labels vector does not match number of samples!\nConsider using the 'subset' argument!\n")
    labels <- newlabels
  } else labels <- data$labels[subset]
  if (!(class %in% names(table(labels))))
    stop(paste(class,"is no class in labels vector!\n"))
  stopifnot(chromosome %in% unique(data$chromosome))

  ### evaluate arguments ###
  geneID <- data$geneName
  geneLoc <- data$geneLocation
  chrom <- data$chromosome
  X <- data$expr[,drop=FALSE]
  if (is.null(kernelparams)){ # set default kernel parameters based on chosen kernel function
    kernelparams <- fitkernelparams(data,chromosome,kernel)
  } #then (is.null(kernelparams))
  isClass <- as.numeric(data$labels == class)
  samResult <- scoring(X,isClass,method="SAM",pcompute="none",verbose=verbose)
  classScores <- matrix(samResult$observed,ncol=1)
  rownames(classScores) <- rownames(X)
  # focus only on one chromosome:
  onChromIndex <- which(chrom==chromosome)
  onChromGenes <- geneID[onChromIndex]
  onChromLocs <- geneLoc[onChromIndex]
  absOnChromLocs <- abs(onChromLocs)
  onChromChrom <- chrom[onChromIndex]
  onChromScores <- classScores[onChromGenes,,drop=FALSE]

  if (is.null(kernelparams))
    kernelparams <- fitkernelparams(data,chromosome,kernel)

  result <- evalParams(locations=absOnChromLocs,scores=onChromScores,kernel=kernel,kernelparams=kernelparams,paramMultipliers=paramMultipliers,ncross=10,verbose=verbose)

  return(result)
} # evaluateParameters

plot.MACATevP <- function(x,...){
  print(plot(x[[1]],x[[2]],type="l",lwd=2,col="red",
             xlab=names(x)[1],ylab="Average Residual Sum-of-Squares",
             main=paste("Cross-validation Results"),...)
        )#print
} # plot.MACATevP
