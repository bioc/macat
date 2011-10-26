\name{discretizeAll}
\alias{discretizeAll}
\title{Discretize complete expression matrix}
\description{
  perform discretization on all chromosomes and write python flat files
}
\usage{
discretizeAll(data, margin = 10)
}
\arguments{
  \item{data}{ MACATData Object }
  \item{margin}{ symmetric quantile in percent }
}
\details{
  The filename for the python flat files are 
  \code{discrete_seqs_margin_<margin>_chrom_<chrom>.py} 
  where <chrom> and <margin>
  are the names of the chromosome and the margin used for discretization.
  For details on the discretization see \code{\link{discretize}}.
}
\value{
  produces python flat file
}
\author{ The MACAT Development team }
\seealso{ \code{\link{discretize}} }
\examples{ 
 # !!! takes some time !!!
 \dontrun{
  #loaddatapkg("stjudem")
  #data(stjude)
  data(stjd)
  discretizeAll(stjd, margin=10) 
 }
}
\keyword{internal}
