\name{preprocessedLoader}
\alias{preprocessedLoader}
\alias{Loader}
\alias{MACATinput}
\alias{MACATData}
\title{Read in data and produce MACAT list}
\description{
 This function reads expression data either from a saved R-file (.RData,.rda),
 or from a tab-separated text-file (.xls). For building a MACAT-list
 structure from objects in your workspace, you can either use this function
 or the convenience wrapper 'buildMACAT'.
}
\usage{
preprocessedLoader(rdatafile, chip, labels = NULL, chromLocObj = NULL,
rdafile = TRUE, tabfile = FALSE, labelfile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdatafile}{Complete name of the expression data file, or the 
	expression matrix}
  \item{chip}{Identifier of the used microarray. To date only commercial
	Affymetrix microarrays are supported by MACAT}
  \item{labels}{Classlabels of the samples, vector of same length as number of
	columns in expression matrix; alternatively complete name of textfile
	with one label per line}
  \item{chromLocObj}{Object of class \code{chromLocation} specifying the
        genomic position, each probe on the array is mapped to. If not
        provided, it is build in the function using \code{annotate}'s
        function \code{buildChromLocation}.}
  \item{rdafile}{Logical; is first argument a saved R-file?}
  \item{tabfile}{Logical; is first argument a tab-separated text file?}
  \item{labelfile}{Logical; is third argument a file with one label per line?}
}
%\details{~~ If necessary, more details than the __description__  above ~~}
\value{
  List of class 'MACATData' with 6 components:
  \item{geneName}{Identifiers of genes/probe sets in expression data}
  \item{geneLocation}{Location of genes on their chromosome as distance from
    5'end in base pairs Negative numbers denote genes on the antisense strand.}
  \item{chromosome}{Chromosome of the respective gene. Components 'geneName',
    'geneLocation', and 'chromosome' are in the same order.}
  \item{expr}{expression matrix with rows = genes and columns = samples/patients}
  \item{labels}{(disease) subtype of each sample, has length = number
	of columns of expression matrix}
  \item{chip}{Identifier for Microarray used for the experiments}
}
%\references{ ~put references to the literature/web site here ~ }
\author{MACAT development team}
\note{At present, \code{macat} can only work with Affymetrix
  microarrays, for which an annotation package is installed on your
  system. Such annotation packages can either be obtained from the
  Bioconductor annotation packages repository or be constructed using
  the Bioconductor package \code{AnnBuilder}. For an example, see the
  common annotation package \code{hgu95av2}.
}
\seealso{\code{\link{buildMACAT}},\code{\link[base]{read.table}},\code{\link{stjd}},\code{stjude} in package 'stjudem'}
\examples{
 \dontrun{
  # assume you have your HG-U95Av2 expression values in a
  # tab-separated text file, called 'foo.txt'
  mydata <- preprocessedLoader("foo.txt","hgu95av2",rdafile=FALSE,tabfile=TRUE)
 }
}
\keyword{file}% at least one, from doc/KEYWORDS
