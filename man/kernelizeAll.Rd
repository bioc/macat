\name{kernelizeAll}
\alias{kernelizeAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth expression data for all chromosomes}
\description{
	'kernelizeAll' smoothes complete expression matrix
	and writes the result into  one
	text file for each chromosome. These text files
	can be read into PYTHON.
}
\usage{
kernelizeAll(data, step.width = 1e+05, kernel = rbf, 
             kernelparams = list(gamma = 1/10^13))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ MACATData Object}
  \item{step.width}{ size of steps for kernelization }
  \item{kernel}{ kernel function one of rbf, kNN, basePairDistance or your owm }
  \item{kernelparams}{ list of named kernel parameters }
}
\details{
  filename of the python flatfiles: 
  \code{kernelized_seq_chrom_<chrom>.py} 
  where <chrom> is the name of the chromosome.
}
\value{
produces text files
}
%\references{ }
\author{ The MACAT Development Team }
%\note{  }

\seealso{ \code{\link{pydata}},  \code{\link{kernelizeToPython}} }
\examples{
 \dontrun{
  # !!! take some time !!!
  loaddatapkg("stjudem")
  kernelizeAll(stjude)
 }
}
\keyword{ manip }% at least one, from doc/KEYWORDS









