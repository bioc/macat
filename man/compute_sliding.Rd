\name{compute.sliding}
\alias{compute.sliding}
\alias{plotSliding}
\title{Compute and plot smoothing of expression values or scores along the chromosome}
\description{
'compute.sliding' computes a smoothing of the expression data or scores
along the chromosome using the specified kernel function. This function
is also used within the 'evalScoring' function.
'plotSliding' creates a plot of the smoothed expression values / scores.
}
\usage{
compute.sliding(data, chromosome, sample, kernel, kernelparams=NULL, step.width = 1e+06)
plotSliding(data, chromosome, sample, kernel, kernelparams=NULL,
step.width=1000000, ...)
}
\arguments{
  \item{data}{ A MACATData list holding the Expression values and gene locations}
  \item{chromosome}{the chromosome to be smoothed}
  \item{sample}{ the sample (patient) whose expression values are smoothed }
  \item{kernel}{ a kernel function (one of rbf, kNN, basePairDistance or your own) }
  \item{kernelparams}{ a list of named parameters for the kernel (by default estimated from the data) }
  \item{step.width}{ the smoothing is computed stepwise every step.width basepairs 
    (default is 100000) }
  \item{...}{further graphical parameters passed on to \code{\link{plot.default}}}
}
\value{
  for \code{compute.sliding}: a matrix of dimension (steps x 2) with
  in the first column the locations in basepairs where an interpolation
  is computed, and in the second column the smoothed values.
  \code{plotSliding} does not return anything and is merely called for
  its side-effect producing the plot.
}
\author{MACAT development team }
\seealso{ \code{\link{kernelize}}, \code{\link{evalScoring}} }
\examples{
data(stjd)
# just compute smoothed values:
smooth = compute.sliding(stjd, chromosome=3, sample=6, rbf,
                         kernelparams=list(gamma=1/10^13))
# compute and plot smoothed values:
plotSliding(stjd, chromosome=3, sample=6,rbf,
            kernelparams=list(gamma=1/10^13),pch=20,
            main="Chromosome 3")
}
\keyword{ manip }
