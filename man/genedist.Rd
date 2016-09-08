\name{genedist}
\alias{genedist}
\title{Chromosomal distribution of the genes in the selected geneset}
\description{The function computes the chromosome wise distribution of the genes in the selected geneset and also plots the chromosomal distribution.}
\usage{genedist(geneset, genelist, plot)}
\arguments{
  \item{geneset}{geneset is a vector of characters representing the names of genes/ gene ids selected from the whole gene list by using a gene selection method.}
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene space: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.}
  \item{plot}{plot is a character string indicating whether the chromosomal distribution of the genes in the selected geneset will be plotted or not. It can be either TRUE/FALSE.}
}
\value{The function returns the chromosomal distribution of the genes in the selected geneset.}
\author{Samarendra Das}
\examples{
data(rice_salt)
data(genelist)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
geneset= GeneSelect(x, y, s=50, method="t-score")$selectgenes
genelist=as.data.frame(genelist)
genedist(geneset, genelist, plot=TRUE)
}
\keyword{gene}
\keyword{geneset}
\keyword{genelist}