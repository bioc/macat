\name{evaluateParameters}
\alias{evaluateParameters}
\alias{plot.MACATevP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluate Performance of Kernel Parameters by Cross-validation}
\description{
  For a given data set, chromosome, class, and kernel function, this function
  helps in determining optimal settings for the kernel parameter(s). The
  performance of individual parameter setting is assessed by cross-
  validation.
}
\usage{
evaluateParameters(data, class, chromosome, kernel, kernelparams = NULL,
                   paramMultipliers = 2^(-4:4), subset = NULL, 
                   newlabels = NULL, ncross = 10, verbose = TRUE)
}
\arguments{
  \item{data}{Gene expression data in the MACAT list format. See data(stjude)
  for an example.}
  \item{class}{Sample class to be analyzed}
  \item{chromosome}{Chromosome to be analyzed}
  \item{kernel}{Choose kernel to smooth scores along the chromose. Available
	are 'kNN' for k-Nearest-Neighbors, 'rbf' for radial-basis-function
	(Gaussian), 'basePairDistance' for a kernel, which averages over
	all genes within a given range of base pairs around a position.}
  \item{kernelparams}{Additional parameters for the kernel as list, e.g., 
	kernelparams=list(k=5) for taking the 5 nearest neighbours in the
	kNN-kernel. If NULL some defaults are set within the function.}
  \item{paramMultipliers}{Numeric vector. If you do cross-validation of the 
	kernel parameters, specify the multipliers of the given (standard)
	parameters to search over for the optimal one.}
  \item{subset}{If a subset of samples is to be used, give vector of column-
	indices of these samples in the original matrix here.}
  \item{newlabels}{If other labels than the ones in the MACAT-list-structure
	are to be used, give them as character vector/factor here. Make sure
	argument 'class' is one of them.}
  \item{ncross}{Integer. Specify how many folds in cross-validation.}
  \item{verbose}{Logical. Should progress be reported to STDOUT?}
}
\value{
  A list of class 'MACATevP' with 4 components:
  \item{[parameterName]}{List of assessed settings for the parameter
                         [parameterName].}
  \item{avgResid}{Average Residual Sum of Suares for the parameter settings in
	the same order as the first component.}
  \item{multiplier}{Multiplier of the original parameters in
	the same order as the first components.}
  \item{best}{List of parameter settings considered optimal by cross-
	validation. Can be directly inserted under the argument
	'kernelparams' of the 'evalScoring' function.}
}
\author{MACAT development team}
\seealso{\code{\link{evalScoring}}}
\examples{
loaddatapkg("stjudem")
evalkNN6 <- evaluateParameters(stjude, class="T", chromosome=6,kernel=kNN, 
                               paramMultipliers=c(0.01,seq(0.2,2.0,0.2),2.5))
plot(evalkNN6)
}
\keyword{manip}% at least one, from doc/KEYWORDS

