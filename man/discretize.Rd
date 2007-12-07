\name{discretize}
\alias{discretize}
\alias{discretizeChromosome}
\alias{discretizeOne}
\title{Discretize expression values}
\description{
  'discretize' returns the discretized expression data for all
  chromosomes in chrom and all samples that have a label listed in label.
  Discretization is performed by comparing the value gene-wise (location-wise) with the symmetric 
  upper and lower quantile given by margin (in percent margin/2 lower and upper quantile).
}
\usage{
discretize(data, chrom, label, margin = 10)
discretizeChromosome(data, chrom, margin=10)
discretizeOne(data, chrom, sample, margin=10)
}
\arguments{
  \item{data}{ MACATData object }
  \item{chrom}{ list of chromosomes }
  \item{label}{ list of labels }
  \item{margin}{ symmetric quantile in percent }
  \item{sample}{ the sample for which you want discretized expression data }
}
\value{
  returns a discretized expression matrix for all genes on the chromosomes
  in 'chrom' and all samples that have a label in 'label'.
}
\author{MACAT development team}
\seealso{\code{\link{discretizeAll}}  }
\keyword{manip}
