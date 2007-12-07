\name{kernelizeToPython}
\alias{kernelizeToPython}
\title{Smooth expression values and write to file}
\description{
  Smoothes expression data for one chromosome and writes 
  the result into a text file, which can be read into PYTHON,
  or returns it without writing.
  'kernelizeToPython' is the one-chromosome version of
  \code{\link{kernelizeAll}}.
}
\usage{
kernelizeToPython(data, chrom, step.width = 1e+05, kernel = rbf, 
                  kernelparams = list(gamma = 1/10^13), saveToFile = TRUE)
}
\arguments{
  \item{data}{ MACATData Object }
  \item{chrom}{ kernelize all genes that are on this chromosome }
  \item{step.width}{ widtrh of interploation steps }
  \item{kernel}{ kernel function one of rbf, kNN, basePairDistance or your own}
  \item{kernelparams}{ list of named kernel parameters }
  \item{saveToFile}{ logical indicating wether to save as flat file or not }
}
\details{
  filename of the flatfile: \code{kernelized_seq_chrom_<chrom>.py}
  where <chrom> is the name of the chromosome.  
}
\value{
  returns kernelized expression matrix
}
\author{ The MACAT Development team }
\seealso{ \code{\link{pydata}}, \code{\link{kernelizeAll}} }
\examples{
  \dontrun{
   data(stjd)
   kernelized = kernelizeToPython(stjd, 3)
  }  
}
\keyword{internal}% at least one, from doc/KEYWORDS
