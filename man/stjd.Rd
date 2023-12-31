% --- Source file: man/stjude.Rd ---
\name{stjd}
\alias{stjd}
\docType{data}
\title{Subset Microarray Data from St.Jude Children Research Hospital (USA)}
\description{
    Example for list-structure used by many functions in MACAT. It's based
    on the gene expression data published by Yeoh et al. (2002)
    The data has been preprocessed using 'vsn' on probe level and the
    probe values have been summed up to probe set values using the
    'median polish' procedure.
    This is a subset of the data, containing only the data for the 5000
    probe sets with the highest variance across the samples and for
    10 exemplary samples, 5 from T-lymphocytic Acute Lymphocytic
    Leukemia (ALL) and 5 from B-lymphocytic ALL.
}
\usage{data(stjd)}
\format{
 List of class 'MACATData' with 6 components:
 \describe{
    \item{geneName:}{Identifiers of genes/probe sets in expression data}
    \item{geneLocation:}{Location of genes on their chromosome as distance from
    5'end in base pairs Negative numbers denote genes on the antisense strand.}
    \item{chromosome:}{Chromosome of the respective gene. Components
      'geneName', 'geneLocation', and 'chromosome' are in the same order.}
    \item{expr:}{expression matrix with rows = genes and columns = samples/patients}
    \item{labels:}{(disease) subtype of each sample, has length = number
      of columns of expression matrix}
    \item{chip:}{Identifier for Microarray used for the experiments
      (here for the Affymetrix HG-U95av2 Oligonucleotide GeneChip)}
  }
}
\note{
  For the full data package see the Bioconductor data package
  \code{stjudem}. If it is not already installed on your system, try
  \code{if (!requireNamespace("BiocManager", quietly=TRUE));
    install.packages("BiocManager");
    BiocManager::install("stjudem")}
}
\references{
    Yeoh et al. Classification, 
    subtype discovery, and prediction of outcome in pediatric acute 
    lymphoblastic leukemia by gene expression profiling. Cancer Cell. 
    March 2002. 1: 133-143.  
}
\seealso{\code{\link{buildMACAT}},\code{stjude} in package 'stjudem' for
the complete expression data}
\examples{
data(stjd)
summary(stjd)
}
\keyword{datasets}
