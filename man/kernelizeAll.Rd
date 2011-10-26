\name{kernelizeAll}
\alias{kernelizeAll}
\title{Smooth expression data for all chromosomes}
\description{
  'kernelizeAll' smoothes complete expression matrix and writes the
  result into  one text file for each chromosome. These text files can
  be read into Python. 
}
\usage{
kernelizeAll(data, step.width = 1e+05, kernel = rbf, 
             kernelparams = list(gamma = 1/10^13))
}
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
  does not return anything; called for its side-effect that is to
  produce Python-readable text files
}
\author{ The MACAT Development Team }
\seealso{ \code{\link{pydata}},  \code{\link{kernelizeToPython}} }
\examples{
 \dontrun{
  # !!! takes quite some time !!!
  loaddatapkg("stjudem")
  data(stjude)
  kernelizeAll(stjude)
 }
}
\keyword{internal}% at least one, from doc/KEYWORDS
