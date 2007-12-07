\name{getResults}
\alias{getResults}
\title{Access results of 'evalScoring'}
\description{
  This function processes the result of the \code{evalScoring} function
  and returns a list of probe sets within chromosome regions deemed
  significant by MACAT. Additional annotation for these probe sets is
  provided along with their identifiers.
}
\usage{
getResults(MACATevalScoringOBJ)
}
\arguments{
  \item{MACATevalScoringOBJ}{Object of class \code{MACATevalScoring},
    usually the result from \code{evalScoring}}
}

\value{
  A list with the following components, describing probe sets within
  chromosome regions deemed significant:
  \item{probeID}{IDs of probe sets within these chromosome regions}
  \item{cytoband}{chromosomal bands these probe sets have been annotated
  to}
\item{geneSYM}{gene symbols these probe sets have been annotated to}
  \item{pvalue}{p-values for probe sets; see details}
  \item{locusid}{EntrezGene-(formerly LocusLink) IDs of these probe sets}
  \item{genedescription}{Description of genes the probe sets have been
    annotated to}
  \item{probeScore}{the differential expression scores for the probe sets}
  \item{chromosome}{chromosome, the analysis has been done for}
  \item{class}{sample class, the analysis has been done for}
}
\author{MACAT development team}
\details{
  The p-values have been computed individually for probe sets (genes),
  not for whole chromosome regions. Thus, regions deemed significant by
  sliding window approach do not have to consist only of probe sets with
  low p-values. These probe-set p-values are not used to determine
  whether a region is considered significant or not. Instead the
  comparison between actual and interpolated scores to actual and
  interpolated boundaries determines whether a region is considered
  significant.

  This function is called within the plot function for the results of
  \code{evalScoring}, when HTML output is desired.
}
\seealso{\code{\link{evalScoring}}}
\examples{
data(stjd)
myevalres <- evalScoring(stjd, class="T", chromosome=6, nperms=10,
                         cross.validate=FALSE) 
results <- getResults(myevalres)
summary(results)
results$probeID[1:20]
}
\keyword{manip}% at least one, from doc/KEYWORDS

