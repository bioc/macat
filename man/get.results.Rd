\name{get.results}
\alias{get.results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Access results of 'evalScoring'}
\description{
  Function returns a list with important information
  about significant chromosome stretches depending on
  a \code{MACATevalScoring} object.
}
\usage{
get.results(MACATevalScoringOBJ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MACATevalScoringOBJ}{MACATevalScoring object}
  }
}
\details{
  Result function for MACATevalScoring object.
}
\value{
  \item{probeID}{Probeset IDs for genes.}
  \item{cytoband}{Cytobands for probeIDs.}
  \item{geneSYM}{GeneSymbol for probeIDs.}
  \item{pvalue}{pValues for probeIDs.}
  \item{locusid}{locusIDs for probeIds.}
  \item{genedescription}{Genedescription for probeIDs.}
  \item{chromosome}{Considered chromosome.}
  \item{class}{Class of the experiments.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{MACAT development team}
%\note{ ~~further notes~~ }

\seealso{\code{\link{evalScoring}}}
\examples{
loaddatapkg("stjudem")
myevalres <- evalScoring(stjude, class="T", chromosome=6, nperms=10,
                         cross.validate=FALSE) 
results <- get.results(myevalres)
summary(results)
results$probeID[1:20]
}

\keyword{manip}% at least one, from doc/KEYWORDS

