\name{evalScoring}
\alias{evalScoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Score differential expression, assess significance,
       and smooth scores along the chromosome}
\description{
  This function computes for all genes on one chromosome the regularized
  t-statistic to score differential gene expression for two given groups
  of samples. Additionally these scores are computed for a number of
  permutations to assess significance. Afterwards these scores are smoothed
  with a given kernel along the chromosome to give scores for chromosomal
  regions.
}
\usage{
evalScoring(data,class,chromosome,nperms=1000,permute="labels",subset=NULL,
     newlabels=NULL,kernel=rbf,kernelparams=NULL,cross.validate=TRUE,
     paramMultipliers=2^(-4:4),ncross=10,step.width=100000,
     memory.limit=TRUE)
}
\arguments{
  \item{data}{Gene expression data in the MACAT list format. See data(stjude)
  for an example.}
  \item{class}{Which of the given class labels is to be analyzed}
  \item{chromosome}{Chromosome to be analyzed}
  \item{nperms}{Number of permutations}
  \item{permute}{Method to do permutations. Default 'labels' does permutations
	of the class labels, which is the common and faster way to assess
	significance of differential expression. The altenative 'locations'
	does permutations of gene locations, is much slower and right only
	preliminary.}
  \item{subset}{If a subset of samples is to be used, give vector of column-
	indices of these samples in the original matrix here.}
  \item{newlabels}{If other labels than the ones in the MACAT-list-structure
	are to be used, give them as character vector/factor here. Make sure
	argument 'class' is one of them.}
  \item{kernel}{Choose kernel to smooth scores along the chromose. Available
	are 'kNN' for k-Nearest-Neighbors, 'rbf' for radial-basis-function
	(Gaussian), 'basePairDistance' for a kernel, which averages over
	all genes within a given range of base pairs around a position.}
  \item{kernelparams}{Additional parameters for the kernel as list, e.g., 
	kernelparams=list(k=5) for taking the 5 nearest neighbours in the
	kNN-kernel. If NULL some defaults are set within the function.}
  \item{cross.validate}{Logical. Should the paramter settings for the kernel
	function be optimized by a cross-validation?}
  \item{paramMultipliers}{Numeric vector. If you do cross-validation of the 
	kernel parameters, specify the multipliers of the given (standard)
	parameters to search over for the optimal one.}
  \item{ncross}{Integer. If you do cross-validation, specify how many folds.}
  \item{step.width}{Defines the resolution of smoothed scores on the
	chromosome, is in fact the distance in base pairs  between 2 
          positions, for which smoothed scores are to be calculated.}
  \item{memory.limit}{If you have a computer with lots of RAM,
	setting this to FALSE will increase speed of computations.}
}

\details{See the documentation for detailed theory and explanations  of this functions.}

\value{
  List of class 'MACATevalScoring' with 11 components:  
  \item{original.geneid}{Gene IDs of the genes on the chosen chromosome, sorted
        according to their position on the chromosome}
  \item{original.loc}{Location of genes on chromosome in base pairs from 5'end}
  \item{original.score}{Regularized t-score of genes on chromosome}
  \item{original.pvalue}{Empirical p-value of genes on chromosome. How often
        was a higher score observed than this one with random permutations?  
        In other words, how significant seems this score to be?}
  \item{steps}{Positions on the chromosome in bp from 5', for which smoothed  
	scores have been computed.}
  \item{sliding.value}{Smoothed regularized t-scores at step-positions.}
  \item{lower.permuted.border}{Smoothed scores from permutations, lower  
	significance border, currently 2.5\%-quantile of permutation scores.}
  \item{upper.permuted.border}{Smoothed scores from permutations, upper  
	significance border, currently 97.5\%-quantile of permutation scores.}
  \item{chromosome}{Chromosome, which has been analyzed}
  \item{class}{Class, which has been analyzed}
  \item{chip}{Identifier for used microarray}
}

\author{MACAT development team}
%\note{}
\seealso{\code{\link{scoring}},\code{\link{plot.MACATevalScoring}},
	\code{\link{getResults}}
}
\examples{
   
    #loaddatapkg("stjudem")
    data(stjd)

    # T-lymphocyte versus B-lymphocyte on chromosome 1, 
    #  smoothed with k-Nearest-Neighbours kernel(k=15), 
    #  few permutations for higher speed
    chrom1Tknn <- evalScoring(stjd,"T",chromosome="1",permute="labels",
    nperms=100,kernel=kNN,kernelparams=list(k=15),step.width=100000)

    # plotting on x11:
    plot(chrom1Tknn)

    # plotting on HTML:
    if (interactive())
       plot(chrom1Tknn,"html")
}
\keyword{manip}% at least one, from doc/KEYWORDS
