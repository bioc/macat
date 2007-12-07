
discretize.tscores <- function(scores) {
  attach(scores)
  sliding.value <- as.vector(sliding.value)
  upper.permuted.border <- as.vector(upper.permuted.border)
  lower.permuted.border <- as.vector(lower.permuted.border)
  detach(scores)
  isUp <- (sliding.value > upper.permuted.border)
  isDown <- (sliding.value < lower.permuted.border)
  isZero <- ((sliding.value < upper.permuted.border) &
             (sliding.value > lower.permuted.border))
  discrete = numeric(length(sliding.value))
  discrete[isUp] = 1
  discrete[isDown] = -1
  discrete[isZero] = 0
  return(discrete)
}

discretizeAllClasses.tscores <- function(data, chrom, nperms=10, kernel=rbf,
                                         kernelparams=NULL, step.width=100000){
  for (label in levels(as.factor(data$labels))) {
    scores = evalScoring(data, label, chromosome=chrom,
      permute="labels",nperms=nperms,subset=NULL,
      newlabels=NULL,kernel=rbf,kernelparams=kernelparams,
      step.width=step.width)
    discrete = discretize.tscores(scores)
    filename = paste("discrete_chrom_", chrom, "_class_", label, ".py",sep="")
    saveForPython(discrete, filename)
  }
}

toPython <- function(list) {
  return(paste("[", paste(list, collapse=","), "]"))
}

discretizeAll <- function(data, margin=10) {
  chromosomes = c(as.character(seq(1, 22)), "X", "Y")
  for(chrom in chromosomes) {
    sink(file=paste("discrete_seqs_margin_", margin,
           "_chrom_", chrom, ".py", sep=""))
    labels = levels(as.factor(data$labels))
    discretized = discretize(data, chrom, data$labels)
    string = list()
    for (label in labels) {
      selected = discretized[, (data$labels == label)]
      pythonMatrix = apply(selected, 2, toPython)
      string = c(string, list(
        paste("[", paste(pythonMatrix, collapse=","), "]\n")))
    }
    cat("[", paste(string, collapse=","), "]\n")
    sink()
  }
}

rawDataToPython <- function(data) {
  chromosomes = c(as.character(seq(1, 22)), "X", "Y")
  for(chrom in chromosomes) {
    cat("Writing chromosome", chrom, "\n")
    sink(file=paste("raw_seqs_chrom_", chrom, ".py", sep=""))
    labels = levels(as.factor(data$labels))
    expr = getExpressionByChromosome.MACATData(data, chrom)
    string = list()
    for (label in labels) {
      selected = expr[, (data$labels == label)]
      pythonMatrix = apply(selected, 2, toPython)
      string = c(string, list(
        paste("[", paste(pythonMatrix, collapse=","), "]\n")))
    }
    cat("[", paste(string, collapse=","), "]\n")
    sink()
  }
}

kernelizeAll <- function(data, step.width=100000, kernel=rbf,
                         kernelparams=list(gamma=1/10^13)){
  for (chrom in c(seq(1, 22), "X", "Y")) {
    cat(paste("Kernelizing chromosome", chrom, "\n"))
    k = kernelizeToPython(data, chrom, step.width, kernel, kernelparams)
  }
}
  
discreteKernelize <- function(data, chrom, margin=10, step.width=100000,
                              kernel=rbf, kernelparams=list(gamma=1/10^13),
                              saveToFile=FALSE) {
  kernelizedExpression = kernelizeToPython(data, chrom, step.width, kernel,
    kernelparams, FALSE)
  upperMargin = apply(kernelizedExpression, 2, quantile, 1 - (margin/200))
  lowerMargin = apply(kernelizedExpression, 2, quantile, margin/200)
  discretized = matrix(1, nrow=dim(kernelizedExpression)[1],
    ncol=dim(kernelizedExpression)[2])
  discretized[kernelizedExpression > upperMargin] = 2
  discretized[kernelizedExpression < lowerMargin] = 0
  if (saveToFile) {
    filename = paste("discrete_kernelized_seq_margin_", margin,
      "_chrom_", chrom, sep="")
    save(discretized, file=paste(filename, ".rdata", sep=""))
    # write to a python file
    cat("Writing to python file\n")
    sink(file=paste(filename, ".py", sep=""))
    string = list()
    for (label in levels(as.factor(data$labels))) {
      selected = discretized[(data$labels == label), ]
      pythonMatrix = apply(selected, 1, toPython)
      string = c(string, list(
        paste("[", paste(pythonMatrix, collapse=","), "]\n")))
    }
    cat("[", paste(string, collapse=","), "]\n")
    sink()
  }
  return(discretized)  
}

kernelizeToPython <- function(data, chrom, step.width=100000, kernel=rbf,
                       kernelparams=list(gamma=1/10^13), saveToFile=TRUE){
  # compute kernelized expressionvalues and discretize them
  genes = (data$chromosome == chrom)
  geneLocations = abs(as.numeric(data$geneLocation[genes]))
  expr = getExpressionByChromosome.MACATData(data, chrom)
  steps = getsteps(geneLocations, step.width)
  cat("Computing kernel weights\n")
  kernelweights = kernelmatrix(steps, geneLocations, kernel, kernelparams)
  cat("Kernelizing\n")
  kernelized = kernelize(expr, kernelweights)
  if (saveToFile == TRUE) {
    filename = paste("kernelized_seq_chrom_", chrom, sep="")
    save(kernelized, file=paste(filename, ".rdata", sep=""))
    # write to a python file
    cat("Writing to python file\n")
    sink(file=paste(filename, ".py", sep=""))
    string = list()
    for (label in levels(as.factor(data$labels))) {
      selected = kernelized[(data$labels == label), ]
      pythonMatrix = apply(selected, 1, toPython)
      string = c(string, list(
        paste("[", paste(pythonMatrix, collapse=","), "]\n")))
    }
    cat("[", paste(string, collapse=","), "]\n")
    sink()
  }
  return(kernelized)
}

discretize <- function(data, chrom, label, margin=10) {
  # this function discretizes the expression values for all sample on one
  # chromosome. discretization takes the expression of each gene over all
  # samples and returns one if the expression is in the upper or lower
  # margin (percent)
  
  all.expr = getExpression.MACATData(data, chrom, label)
  upperMargin = apply(all.expr, 1, quantile, 1 - (margin/200))
  lowerMargin = apply(all.expr, 1, quantile, margin/200)
  discretized = matrix(1, nrow=dim(all.expr)[1], ncol=dim(all.expr)[2])
  discretized[all.expr > upperMargin] = 2
  discretized[all.expr < lowerMargin] = 0
  return(discretized)  
}
  
discretizeChromosome <- function(data, chrom, margin=10) {
  # this function discretizes the expression values for all sample on one
  # chromosome. discretization takes the expression of each gene over all
  # samples and returns one if the expression is in the upper or lower
  # margin (percent)
  
  all.expr = getExpressionByChromosome.MACATData(data, chrom)
  upperMargin = apply(all.expr, 1, quantile, 1 - (margin/200))
  lowerMargin = apply(all.expr, 1, quantile, margin/200)
  discretized = matrix(0, nrow=dim(all.expr)[1], ncol=dim(all.expr)[2])
  discretized[all.expr > upperMargin] = 1
  discretized[all.expr < lowerMargin] = -1
  return(discretized)  
}

discretizeOne <- function(data, chrom, sample, margin=10) {
  # this function discretizes the expression values for one sample on one
  # chromosome. discretization takes the expression of each gene over all
  # samples and returns one if the expression is in the upper or lower
  # margin (percent)
  
  all.expr = getExpressionByChromosome.MACATData(data, chrom)
  this.expr = all.expr[, sample]
  upperMargin = apply(all.expr, 1, quantile, 1 - (margin/200))
  lowerMargin = apply(all.expr, 1, quantile, margin/200)
  discretized = rep(0, length(this.expr))
  discretized[this.expr > upperMargin] = 1
  discretized[this.expr < lowerMargin] = -1
  return(discretized)
}
  

saveForPython <- function (list, filename) {
  sink(file=filename)
  cat(toPython(list))
  sink()
}

