\name{kernelize}
\alias{kernelmatrix}
\alias{getsteps}
\alias{kernelize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth expression values or scores}
\description{
  'kernelize' uses a kernel to smooth the data given in geneLocations by
  computing a weighted sum of the values vector. The weights for each position
  are given in the kernelweights matrix. A kernelweights matrix can be obtained
  by using the kernelmatrix function.
}
\usage{
getsteps(geneLocations, step.width)
kernelmatrix(steps, geneLocations, kernel, kernelparams)
kernelize(values, kernelweights)
}
\arguments{
	\item{geneLocations}{ a list of gene locations (length n)}
	\item{step.width}{ the width of steps in basepairs }
	
	\item{steps}{ a list of locations where the kernelization shall be computed }
	\item{kernel}{kernel function one of rbf, kNN or basePairDistance (or your own) }
	\item{kernelparams}{ a list of named parameters for the kernel (default is fitted to the data) }
	
	\item{values}{vector of length n or matrix (m x n) of values that are to be smoothed}
	\item{kernelweights}{a matrix of (n x steps) where n is the length of the
  		values vector and steps is the number of points where you wish to interpolate }

}
\value{
  \item{getsteps}{a list of locations starting at min(genLocations) going to max(geneLocations) with steps of size step.width}
  \item{kernelmatrix}{ a matrix of (n x steps) containing the kernel weights for each location in steps}
  \item{kernelize}{a vector of length steps or a matrix (m x steps) containing the smoothed values}
}

%\references{ }
\author{ MACAT Development team }
%\note{  }


\seealso{ \code{\link{compute.sliding}}, \code{\link{evalScoring}} }
\examples{
  loaddatapkg("stjudem")
  genes = seq(100)
  geneLocations = abs(stjude$geneLocation[genes])
  geneExpression = stjude$expr[genes,]
  step.width = 100000
  steps = getsteps(geneLocations, step.width)
  weights = kernelmatrix(steps, geneLocations, rbf, list(gamma=1/10^13))
  kernelized = kernelize(geneExpression, weights)
  plot(steps, kernelized[1,])
}
\keyword{ manip }% at least one, from doc/KEYWORDS
