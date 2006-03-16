###############################################################################
# This file holds all functions needed for loading and preprocessing of data
# definition of a general list for chromosomal location and
# expression levels
###############################################################################
###############################################################################
# function preprocessedLoader:
###############################################################################
# loads the preprocessed datamatrix from 'rdatafile' using the annotate package
# to determine the chromosomal locations of the genes and returns a data frame
# the leading columns contain information about the gene and its location, the
# trailing columns hold the expressionlevels for all samples
###############################################################################

preprocessedLoader <- function(rdatafile,chip,labels=NULL,rdafile=TRUE,tabfile=FALSE,labelfile=FALSE) {
  # file: is first argument  a file or the expression matrix
  require(annotate)
  #acceptedChips <-  c("hu6800","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e",
  #                    "hgu133a","hgu133b","mgu74av2","mgu74bv2","mgu74c",
  #                    "moe430a","moe430b","rae230a","rae230b","rgu34a",
  #                    "rgu34b","rgu34c")
  loadChip <- require(chip, character.only=TRUE, quiet=TRUE)
  if (!loadChip){
    stop(paste("Unknown Chip! At present, MACAT can only deal with chips, for which there",
               "exists an installed BioConductor annotation data package, such as 'hgu95av2'.",
               "The chip type must be specified using the name of this annotation package,",
               "e.g., for the Human Genome 95 A version 2 chip: 'hgu95av2'.\n", sep="\n"))
  }
  cat("Reading expression matrix...\n")
  if ((!rdafile)&(!tabfile)) # is m already a matrix?
    m <- as.matrix(rdatafile)
  else if (rdafile) {
    oldworkspace <- c()
    oldworkspace <- ls()
    load(rdatafile)
    newworkspace <- ls()
    newobject <- setdiff(newworkspace,oldworkspace)
    if (length(newobject)!=1)
      stop ("R-Data file contained no single new object!\nTry specifying expression matrix as first argument and set argument 'file=FALSE'!\n")
    m <- eval(as.symbol(newobject))
    #m=expr.mat
  } else if (tabfile){
    m <- as.matrix(read.delim(rdatafile,row.names=1))
  }
  if (!is.null(labels)){
    cat("Reading labels...\n")   
    if (labelfile)
      labels <- as.vector(as.matrix(read.table(labels)))
    if (length(labels)!=ncol(m))
      stop(paste("Third argument 'labels' indicates",length(labels)," samples, while expresssion matrix contains",ncol(m),"samples (=columns)!\n"))
  } # if (!is.null(labels))
                    
  cat("Assessing chromosome information for genes on chip...\n")
  chromLocationObj <- buildChromLocation(chip)
  if (is.null(colnames(m))){
    colnames(m) <- paste("sample",as.character(1:ncol(m)),sep="")
  }
  # build pheno data to satisfy current exprSet validation:
  sampleDataFrame <- data.frame(sampleNames=I(colnames(m)),
                                row.names=colnames(m))
  pdata <- new("phenoData", pData=sampleDataFrame,
               varLabels=list(names(sampleDataFrame)))
  eset  <- new("exprSet", exprs=m, phenoData=pdata)
  chromosomes <- names(chromLocationObj@chromInfo)                  
  allGeneNames <- c()
  allGeneLocations <- c()
  allChromosomes <- c()

  # auxiliary function from 'annotate' with fixed bug:
  usedChromGenes2 <- function (eSet, chrom, specChrom){
    cLocs <- chromLocs(specChrom)
    genes <- cLocs[[chrom]]
    usedGenes <- genes[names(genes) %in% geneNames(eSet)]
    if (length(usedGenes)==0) return(NULL)
    ord <- order(abs(usedGenes))
    usedGenes <- as.list(usedGenes[ord])
    return(usedGenes)
  }#usedChromGenes2
  
  for (chromosome in chromosomes) {
    cat(paste("Locating genes on chromosome",chromosome,".... "))
    usedGenes <- usedChromGenes2(eset, chromosome, chromLocationObj)
    if (length(usedGenes)==0){
      cat("0\n")
      next
    }
    usedGenes <- usedGenes[!is.na(usedGenes)]
    cat(length(usedGenes),"\n")
    if (length(usedGenes)==0) next
    geneName <- names(usedGenes)
    geneLocation <- as.numeric(usedGenes)
    allGeneNames <- c(allGeneNames,geneName)
    allGeneLocations <- c(allGeneLocations,geneLocation)
    allChromosomes <- c(allChromosomes,as.character(rep(chromosome,length(geneName))))
  } # for chromosome
  if (length(allGeneNames)==0)
    stop("Gene-Identifier not found on specified chip!\nIf Affymetrix-Chip,Probe-Set-IDs are expected!")
  data <- list(geneName=allGeneNames,geneLocation=allGeneLocations,
               chromosome=allChromosomes, expr=m,labels=labels,chip=chip)
  class(data) <- "MACATData"
  return(data)  
} # preprocessedLoader

## wrapper for preprocessedLoader: ###
buildMACAT <- function(matrix,chip,labels=NULL){
  if (class(matrix)=="exprSet")
    matrix <- exprs(matrix)
  if (!is.matrix(matrix)) stop("First argument is no matrix!\n")
  return(preprocessedLoader(matrix,chip,labels,rdafile=FALSE,
                            tabfile=FALSE,labelfile=FALSE)) 
} #buildMACAT


###########################################################################
## class method for accessing all expression values of one chromosome
###########################################################################

getExpressionByChromosome.MACATData <- function(data, chrom) {
  return(getExpression.MACATData(data, chrom, data$labels))
}

getExpressionByClass.MACATData <- function(data, label) {
  return(getExpression.MACATData(data, data$chromosome, label))
}

getExpression.MACATData <- function(data, chrom, label) {
  genesOnChrom = (data$chromosome == chrom)
  samples = (data$labels == label)
  geneNamesOnChrom = data$geneName[genesOnChrom]
  expressionsOnChrom = data$expr[geneNamesOnChrom, samples]
  return(expressionsOnChrom)
}

getExpressionByProbeset.MACATData <- function(data, probeset) {
  return(0)
}
