\name{compute.sliding}
\alias{compute.sliding}
\alias{plotSliding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute and plot smoothing of expression values or scores along the chromosome}
\description{
'compute.sliding' computes a smoothing of the expression data or scores
along the chromosome using the specified kernel function. This function
is also used within the 'evalScoring' function.
'plotSliding' creates a plot of the smoothed expression values / scores.
}
\usage{
compute.sliding(data, chromosome, sample, kernel, kernelparams=NULL, step.width = 1e+06)
plotSliding(data, chromosome, sample, kernel, kernelparams=NULL, step.width=1000000)
}
\arguments{
  \item{data}{ A MACATData list holding the Expression values and gene locations}
  \item{chromosome}{the chromosome to be smoothed}
  \item{sample}{ the sample (patient) whose expression values are smoothed }
  \item{kernel}{ a kernel function (one of rbf, kNN, basePairDistance or your own) }
  \item{kernelparams}{ a list of named parameters for the kernel (by default estimated from the data) }
  \item{step.width}{ the smoothing is computed stepwise every step.width basepairs 
  (default is 100000) }
}
\details{
}
\value{
  a matrix of dimension (steps x 2) with
  in the first column the locations in basepairs where an
   interpolation is computed, and
  in the second column the smoothed values.
}
%\references{ }
\author{MACAT development team }
%\note{  }

\seealso{ \code{\link{kernelize}}, \code{\link{evalScoring}} }
\examples{
  loaddatapkg("stjudem")
  smooth = compute.sliding(stjude, chromosome=3, sample=6, rbf, kernelparams=list(gamma=1/10^13))
  plotSliding(stjude, chromosome=3, sample=6, rbf, kernelparams=list(gamma=1/10^13))
}
\keyword{ manip }
