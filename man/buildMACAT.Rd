\name{buildMACAT}
\alias{buildMACAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create MACAT list from objects in workspace}
\description{
	This is a wrapper around the preprocessedLoader function. Use it, 
	when you want to build a MACAT-list structure from objects already
	in your workspace.
}
\usage{
buildMACAT(matrix, chip, labels = NULL, chromLocObj = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{expression matrix with rows=genes and columns=samples;
	Rownames have to match chip; Columnnames are not mandatory.}
  \item{chip}{Identifier for used microarray}
  \item{labels}{Classlabels for samples, has to have length=number of columns
    in matrix}
  \item{chromLocObj}{Object of class \code{chromLocation} specifying the
  genomic position, each probe on the array is mapped to. If not
  provided, it is build in the function using \code{annotate}'s
  function \code{buildChromLocation}.}
}
\details{
  This is only a convenience wrapper around the function
  \code{\link{preprocessedLoader}} for the case, that you want to build
  a MACAT-list from objects in your workspace.
}
\value{
  A MACAT-list structure. For an example and a description of the format
  see data \code{stjude} in package 'stjudem'.
}
\author{MACAT development team}
%\note{ ~~further notes~~ }
\seealso{\code{\link{preprocessedLoader}},\code{stjude} in package 'stjudem'}
\examples{
 X <- matrix(rnorm(200),nrow=20,ncol=10)
 rownames(X) <- c('34916_s_at','34917_at','34462_at','163_at','35219_at',
	'31641_s_at','33300_at','33301_g_at','38950_r_at','41249_at',
	'294_s_at','32004_s_at','33299_at','41243_at','33341_at','362_at',
	'1918_at','41499_at','41500_at','41282_s_at')
 colnames(X) <- paste("Sample",1:10,sep="")
 y <- rep(c("A","B"),c(5,5))
 toy <- buildMACAT(X,"hgu95av2",y)
 summary(toy)
}
\keyword{file}% at least one, from doc/KEYWORDS
