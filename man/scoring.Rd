\name{scoring}
\alias{scoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute (regularized) t-scores for gene expression data}
\description{
 This function computes for all genes in an expression matrix the (regularized)
 t-scores (statistics) with the given class labels and a number of 
 permutations of these labels. Each gene is also assigned a p-value either 
 empirically from the permutation scores or from a t-distribution.
}
\usage{
scoring(data, labels, method = "SAM", pcompute = "tdist", 
        nperms = 1000, memory.limit = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Expression matrix with rows = genes and columns = samples}
  \item{labels}{Vector or factor of class labels; Scoring works only with
	two classes!}
  \item{method}{Either "SAM" to compute regularized t-scores, or "t.test" to
 	compute Student's t-statistic}
  \item{pcompute}{Method to compute p-values for each genes, either
	"empirical" to do permutations and compute p-values from them,
	or "tdist" to compute p-values based on respective t-distribution}
  \item{nperms}{Number of permutations of the labels to be investigated, if
	argument 'pcompute="empirical"'}
  \item{memory.limit}{Logical, if you have a really good computer (>2GB RAM), 
	setting this FALSE will increase speed of computations}
  \item{verbose}{Logical, if progress should be reported to STDOUT}
}
\details{
  If 'pcompute="empirical"', the statistic is computed based on the given class
  labels, afterwards for 'nperms' permutations of the labels. The p-value for
  each gene is then the proportion of permutation statistics that are higher or
  equal than the statistic from the real labels. For each gene the 2.5\%- and
  the 97.5\%-quantile of the permutation statistics are also returned as lower
  and upper 'significance threshold'.

  If 'pcompute="tdist", the statistic is computed only based on the given
  class labels, and the p-value is computed from the t-distribution with
  (Number of samples - 2) degrees of freedom.  
}
\value{
  A list, with four components:
  \item{observed}{(Regularized) t-scores for all genes based on the given labels}
  \item{pvalues}{P-values for all genes, either from permutations or t-distribution}
  \item{expected.lower}{2.5\%-quantile of permutation test-statistics, supposed
	to be a lower 'significance border' for the gene; or NULL if
	p-values were computed from t-distribution}
  \item{expected.upper}{2.5\%-quantile of permutation test-statistics, supposed
	to be an upper 'significance border' for the gene; or NULL if
	p-values were computed from t-distribution}
}
\references{Regarding the regularized t-score see the \code{macat_docu.pdf}}
\author{MACAT development team}
\note{In MACAT this function is only called internally by \code{\link{evalScoring}}}
\seealso{\code{\link{evalScoring}}}
\examples{
 \dontrun{
  data(stjd)
  # compute scores for T- vs. B-lymphocyte ALL:
  isT <- as.numeric(stjd$labels=="T")
  TvsB <- scoring(stjd$expr,isT,method="SAM",pcompute="empirical",nperms=100)
  summary(TvsB$observed)
 }
}
\keyword{manip}% at least one, from doc/KEYWORDS
