\name{discretize.tscores}
\alias{discretize.tscores}
\alias{discretizeAllClasses.tscores}
\title{Discretize regularized t-scores}
\description{
  discretize.tscores returns a discretized version of the scores in the MACATevalScoring object.
  Discretization is performed by comparing the value gene-wise (location-wise) with the symmetric 
  upper and lower quantile given by margin (in percent margin/2 lower and upper quantile).
  discretizeAllClasses produces a flatfile readable by PYTHON.
}
\usage{
discretize.tscores(scores)
discretizeAllClasses.tscores(data, chrom, nperms=10, kernel=rbf, kernelparams=NULL, step.width=100000)
}
\arguments{
  \item{scores}{ a MACATevalScoring object obtained from evalScoring }
  \item{data}{ a MACATData Object containing all expression values, geneLocations and labels 
  (obtained from preprocessedLoader) }
  \item{chrom}{ chromosome that is discretized }
  \item{nperms}{ number of permutations for the computation of empirical p values (evalScoring) }
  \item{kernel}{ kernel function used for smoothing one of rbf, kNN, basePairDistance or your own }
  \item{kernelparams}{ list of parameters for the kernels }
  \item{step.width}{ size of a interpolation step in basepairs }
}
\details{
  The filename for the python flat files are 
  \code{discrete_chrom_<chrom>_class_<label>.py} 
  where <chrom> and <label>
  are the names of the chromosome and class label.
}
\value{
	\item{discretize.tscores}{ a vector of discretized tscores }
	\item{discretizeAllClasses.tscores}{ creates python flatfiles (see details) }
}
\author{  The MACAT development team }
\seealso{ \code{\link{evalScoring}}, \code{\link{kernels}}, \code{\link{pythondata}}}
\examples{
  #loaddatapkg("stjudem")
  #data(stjude)
  data(stjd)
  # simple scoring with short running time
  scores = evalScoring(stjd, "T", 1, nperms=100, cross.validate=FALSE)
  discrete = discretize.tscores(scores)
}
\keyword{manip}
