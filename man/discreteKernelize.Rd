\name{discreteKernelize}
\alias{discreteKernelize}
\title{Discretize and smooth expression values}
\description{
  returns discretized kernelized expression values and saves them to a
  file if argument 'saveToFile' is TRUE. 
  For details on discretization see \code{\link{discretize}}.
}
\usage{
discreteKernelize(data, chrom, margin = 10, step.width = 1e+05, kernel = rbf,
                  kernelparams = list(gamma = 1/10^13), saveToFile = FALSE)
}
\arguments{
  \item{data}{ MACATData Object }
  \item{chrom}{ chromosome to kernelize}
  \item{margin}{ symmetric qunatile in percent }
  \item{step.width}{ size of the interpolation steps }
  \item{kernel}{ kernel function one of rbf, kNN, basePairDistance or your own }
  \item{kernelparams}{ list of named kernel parameters }
  \item{saveToFile}{logicval indicating whether to write a flatfile or
    not; default is FALSE}
}
\details{
  Filename of the flatfile is: 
  \code{discrete_kernelized_seq_margin_<margin>_chrom_<chrom>.py}
  where <margin> is the discretization parameter and <chrom> the name of the
  chromosome.
}
\value{
  discretized and kernelized expression matrix
}
\author{ The MACAT Development team }
\seealso{\code{\link{pydata}}, \code{\link{kernelizeAll}} }
\examples{
   #loaddatapkg("stjudem")
   #data(stjude)
   data(stjd)
   discretizedKernelized = discreteKernelize(stjd, 13)
}
\keyword{manip}
