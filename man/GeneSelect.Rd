\name{GeneSelect}
\alias{GeneSelect}
\title{Selection of informative geneset}
\description{The function returns the informative geneset from the high dimensional gene expression data using a proper statistical technique.}
\usage{GeneSelect(x, y, s, method)}
\arguments{
  \item{x}{x is a N x m gene expression data matrix (must be data frame) and row names as gene names, where, N represents the number of genes in the whole gene space and m is number of samples.}
  \item{y}{y is a m by 1 vector representing the sample labels, is according to the different stress conditions for two class problem (must be 1: stress/-1: control)}
  \item{s}{s is a numeric constant representing the number of genes to be selected from the large pool of genes/ gene space.}
  \item{method}{method is a character string indicating which method for informative gene selection is to be used. One of method "t-score" (default), "F-score", "MRMR", "BootMRMR"  can be abbreviated and used.}
}
\value{The function returns the informative geneset using a particular method from the high dimensional gene expression data.}

\author{Samarendra Das}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
GeneSelect(x, y, s=50, method="t-score")$selectgenes
}
\keyword{gene}
\keyword{geneset}
\keyword{gene expression data}
\keyword{method}
