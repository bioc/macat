\name{plot.MACATevalScoring}
\alias{plot.MACATevalScoring}
\title{Plot function for MACATevalScoring objects.}
\description{
  Function plots scores, 0.025 and 0.975 quantiles of the permuted
  scores (grey lines), and sliding average score (red line) along the
  chromosome. Yellow dots highlight regions, in which the smoothed
  absolute scores exceed the permutation-derived quantile boundaries.
}
\usage{
\method{plot}{MACATevalScoring}(x, output = "x11", 
              HTMLfilename = NULL, mytitle = NULL, new.device = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{MACATevalScoring object.}
  \item{output}{plot "x11" or create a "html" -file with further information. HTML-page will open automatically.}
  \item{HTMLfilename}{HTML-filename, default:\code{Results<CHOMOSOME>_<CLASS>.html}}
  \item{mytitle}{Title of HTML-page, default: \dQuote{Results of class <CLASS> on chromosome <CHROMOSOME>}}
  \item{new.device}{if FALSE: Possibility to plot several plots in one device}
  \item{...}{further arguments passed on to generic function \code{plot}}
}
\details{
  One can create a HTML-page on-the-fly if argument output='html'. 
  The HTML-page provides informations about highlighted regions in the plot. 
  Furthermore there are click-able Entrezgene-IDs for further analysis. 
}
\author{MACAT development team}
\seealso{\code{\link{evalScoring}}, \code{\link{getResults}}}
\examples{
 # see function 'evalScoring' for an example
}
\keyword{hplot}% at least one, from doc/KEYWORDS
