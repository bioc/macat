\name{Kernels}
\alias{kNN}
\alias{rbf}
\alias{basePairDistance}
\alias{fitkernelparams}
\alias{kernels}
\title{various kernel functions for computations in MACAT}
\description{
    Various Kernel functions for computations in MACAT. Normally not
    called by user.
    All kernel functions have the same arguments in the same order!!!
}
\usage{
kNN(geneLocations, position, params)
rbf(geneLocations, position, params = list(gamma=1/10^13))
basePairDistance(geneLocations, position, params = list(distance = 1e+06))
}
\arguments{
  \item{geneLocations}{Location of genes}
  \item{position}{Position on chromosome}
  \item{params}{special kernel parameters
  	\describe{
		\item{kNN:}{k = number of nearest genes}
		\item{rbf:}{gamma = kernel width}
		\item{basePairDistance:}{distance = distance within which all
		genes are averaged }
	}
      }
}
\details{
For internal use by other MACAT-functions. Not called by user.
}
\value{
  returns kernel weight for position, computed from the geneLocations
}
\author{MACAT development team}
\seealso{\code{\link{evalScoring}, \link{compute.sliding}, \link{kernelize}}}
\examples{
  data(stjd)
  genes = seq(100)
  geneLocations = abs(stjd$geneLocation[genes])
  position = c(1000) # location for which you want the kernelweights
  kernelweights = rbf(as.matrix(geneLocations), as.matrix(position), 
  		      list(gamma=1/10e13))
  hist(kernelweights)  
}

\keyword{internal}
