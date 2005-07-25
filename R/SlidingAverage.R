
squaredDist <- function(x,y){
  # x,y: (column) vectors
  z <- as.vector(x-y)
  return(as.vector(t(z) %*% z))
} #squaredDist

#--------------------------------------------------------------------------
rbf <- function(geneLocations, position, params=list(gamma=1/10^13)) {
  # radial basis function kernel for x, y
  # params is a named list: gamma is the width of the kernel
  # default value of gamma (variance) is chosen such that the standard
  # deviation is one MB (->1 cmorgan)
  # argumentation with mean distance between genes could also explain.
  gamma = params$gamma
  kernel<-function(x,y,gamma){return(exp(-gamma * squaredDist(x,y)))}
  return(apply(geneLocations, 1, kernel, position, gamma))
} #rbf


#--------------------------------------------------------------------------
kNN <- function(geneLocations, position, params) {
  # kernel function for k nearest neighbours
  # returns 1 if y is one of the k nearest neighbours of x, 0 otherwise
  # params is a named list: data gives a matrix of all data in rows
  # k gives the number of neirest neighbours
  k = params$k
  #data = params$data
  #data = as.matrix(data)
  dists <- abs(geneLocations-position)
  knn <- sort(dists,method="quick",index.return=TRUE)$ix[1:k]
  weights <- numeric(length(geneLocations))  
  weights[knn] = 1
  return(weights)
} #kNN

#--------------------------------------------------------------------------
basePairDistance <- function(geneLocations, position, params = list(distance = 1000000)) {
  # gene is a gene coordinate in basepairs
  # x: position around which you want to find all neighbours within distance
  # returns 1 if gene is within distance from x
  distance = params$distance
  kernel <- function(geneLocation, position, distance){
    if (abs(geneLocation - position) < distance) {
      return(1)
    }
    return(0)
  }
  return(apply(geneLocations, 1, kernel, position, distance))
}#basePairDistance

#--------------------------------------------------------------------------
plotSliding <- function(data, chromosome, sample, kernel, kernelparams=NULL, step.width=1000000, ...) {
  if (is.null(kernelparams))
    kernelparams=fitkernelparams(data,chromosome,kernel)
  # find expressions of the used genes on the chromosome for the sample
  genesOnChromIndex <- which(data$chromosome == chromosome)
  genesOnChrom <- data$geneName[genesOnChromIndex]
  # take the abs because we dont care on which strand the genes lie
  genes <- abs(data$geneLocation[genesOnChromIndex])
  expr <- data$expr[genesOnChrom,sample]
  points = compute.sliding(data, chromosome, sample, kernel, kernelparams, step.width)
  steps = points[,1]
  sliding.value = points[,2]
  print(plot(genes, expr, "p", ylab="Expression", xlab="Coordinate", ylim=c(min(c(expr, sliding.value,0)), max(expr, sliding.value)),...))
  lines(steps, sliding.value, col="red",lwd=2)
} # plotSliding

#--------------------------------------------------------------------------
compute.sliding <- function(data, chromosome, sample, kernel, kernelparams=NULL, step.width=1000000) {
  if (is.null(kernelparams))
    kernelparams=fitkernelparams(data,chromosome,kernel)
  # find expressions of the used genes on the chromosome for the sample
  genesOnChromIndex <- which(data$chromosome == chromosome)
  genesOnChrom <- data$geneName[genesOnChromIndex]
  
  # take the abs because we dont care on which strand the genes lie
  genes <- abs(data$geneLocation[genesOnChromIndex])
  ngenes <- length(genes)
  expr <- data$expr[genesOnChrom,sample]
  #kernelparams = c(kernelparams, list(data=genes))
  # the coordinates of the used genes are given in used genes
  # slide from the min to the max
  steps <- getsteps(genes,step.width)
  #steps = seq(min(as.numeric(genes)), max(as.numeric(genes)), step.width)
  # compute the slidingAverageAt the positions in step using the usedGenes
  # coordinates and the expression values at these coordinates
  # with the knn kernel and the paramterlist
  kernelweights <- kernelmatrix(steps,genes,kernel,kernelparams)
  sliding.value = t(expr) %*% kernelweights
  return(cbind(steps, t(sliding.value)))
} 
#-----------------------------------------------------------------------
getsteps <- function(geneLocations,step.width=1000000){
  steps <- seq(min(as.numeric(geneLocations)),max(as.numeric(geneLocations)),step.width)
  return(steps)
} #getsteps

#--------------------------------------------------------------------
kernelmatrix <- function(steps,geneLocations,kernel,kernelparams){
  kernelweights <- apply(as.matrix(steps),1,
                         function(x){k = kernel(as.matrix(geneLocations),x,kernelparams);
                                     if (sum(k)==0)
                                       return(numeric(length(geneLocations)))
                                     else
                                       return(k/sum(k))}) 
  return(kernelweights)
} #kernelmatrix
#--------------------------------------------------------------------------
kernelize <- function(values,kernelweights){
  sliding.value <- t(values) %*% kernelweights
} #kernelize

#--------------------------------------------------------------------------
maxDistances <- function(m, chromosomes) {
  distances = numeric()
  for (chrom in chromosomes) {
    #locations = abs(m[m$chromosome == chrom,]$geneLocation) #alter data frame
    genesOnChromIndex <- which(m$chromosome == chrom)
    locations = abs(m$geneLocation[genesOnChromIndex])
    locations.u = unique(locations) 
    locations.u = sort(locations)
    
    for (i in seq(length(locations.u) - 1)) {
      distances = c(distances, locations.u[i+1] - locations.u[i])
    }
  }
  #hist(distances, breaks=100)
  return(max(distances))
}

#--------------------------------------------------------------------------
pair.distances <- function(m, chromosomes) {
  distances = numeric()
  for (chrom in chromosomes) {
    locations = abs(m[m$chromosome == chrom,]$geneLocation) #alter data frame
    locations.u = unique(locations)
    locations.u = sort(locations)
    
    for (i in seq(length(locations.u))) {
      for (j in seq(i, length(locations.u))) {
        distances = c(distances, abs(locations.u[i] - locations.u[j]))
      }
    }
  }
  return(distances)
}


#--------------------------------------------------------------------------
# this function visualizes the effect of different choices of gamma for rbf
# the 
compare.gammas <- function(m, chromosome, sample) {
  max = maxDistances(m, chromosome)
  if (interactive() && capabilities()["X11"])
    x11(width=16, height=8)  
  plotSliding(m, 6, 3, rbf, list(gamma=log(2)/((max/2)^2)))
  maxHalf = compute.sliding(m, 6, 3, rbf, list(gamma=1/((max/2)^2)))
  lines(maxHalf, col="blue", lwd=2)
  default =  compute.sliding(m, 6, 3, rbf, list(gamma=1/(10000000000000)))
  lines(default, col="green", lwd=2)
  legend(0,5.7, c("red log(2)/((max/2)^2)","blue 1/((max^2)/2)", "green 10^13"))
}


# ----------------------------------
fitkernelparams <- function(data,chromosome,kernel){
  if(identical(kernel,rbf)){
    maxDist <- maxDistances(data,chromosome)
    gamma <- log(2)/(maxDist/2)^2
    # see macat_docu.pdf for detailed explanation
    kernelparams <- list(gamma=gamma)
  } # if(identical(kernel,rbf))
  if(identical(kernel,kNN)){
    nGenesOnChrom <- sum(data$chromosome==chromosome)
    kernelparams <- list(k=ceiling(nGenesOnChrom/10))
  }
  if(identical(kernel,basePairDistance)){
    maxDist <- maxDistances(data,chromosome)
    kernelparams <- list(distance=maxDist)
  }
  return(kernelparams)
} # fitkernelparams

