\name{loaddatapkg}
\alias{loaddatapkg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load data package}
\description{
 This function loads the data package, you need for seeing the demo and
 the examples. If you have already installed the data package, it will
 simply attach it via \code{library}. Otherwise it will try to download
 and install the package using functions from the package
 \code{reposTools}. If you have not installed \code{reposTools} either,
 it will prompt you to proceed to a web page, download the package from
 there, and install it manually.
}
\usage{
 loaddatapkg(mydatapkg,installDir=.libPaths()[1])
}
\arguments{
  \item{mydatapkg}{Name of the data package to load as String}
  \item{installDir}{Directory, into which the new package will be
		installed, if is not already installed. Defaults
		to the first entry of the standard installation
		paths.}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
%\value{
%  ~Describe the value returned
%}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling}
%\note{ ~~further notes~~ }
\seealso{\code{\link[base]{library}}, package \code{reposTools}}
\examples{
 \dontrun{loaddatapkg("stjudem") # to load the data package "stjudem"}
}
\keyword{internal}% at least one, from doc/KEYWORDS
