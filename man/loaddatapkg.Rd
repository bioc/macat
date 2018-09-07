\name{loaddatapkg}
\alias{loaddatapkg}
\title{Load data package}
\description{
 This function loads the data package, you need for seeing the demo and
 the examples. If you have already installed the data package, it will
 simply attach it via \code{library}. Otherwise it will try to download
 and install the package using functions from the package
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
\note{The package \code{stjudem} by now is a Bioconductor example data
  package, too. Thus, you can also install with the following:
  \code{if (!requireNamespace("BiocManager", quietly=TRUE));BiocManager::install("stjudem")}
    install.packages("BiocManager");BiocManager::install("stjudem")}
}
\author{Joern Toedling}
\seealso{\code{\link[utils:update.packages]{install.packages}}}
\examples{
 \dontrun{loaddatapkg("stjudem") # to load the data package "stjudem"}
}
\keyword{internal}% at least one, from doc/KEYWORDS
