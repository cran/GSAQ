\name{qtldist}
\alias{qtldist}
\title{QTL wise distribution of genes in the selected geneset}
\description{Computation of number of qtl-hit genes in each QTL and also QTL wise distribution  of genes in the selected geneset }
\usage{qtldist(geneset, genelist, qtl, plot)}
\arguments{
  \item{geneset}{geneset is a vector of characters representing the names of genes/ gene ids selected from the whole gene space by using a gene selection method.}
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene space: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.}
  \item{qtl}{qtl is a Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.}
  \item{plot}{plot is a character string used to plot the QTL wise distribution of genes in the selected gene set. It can be either TRUE/FALSE.}
}
\value{The function returns number of qtl-hit genes in each QTL and QTL wise distribution of the selected genes.}

\author{Samarendra Das}
\examples{
data(rice_salt)
data(genelist)
data(qtl_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
geneset= GeneSelect(x, y, s=50, method="t-score")$selectgenes
genelist=as.data.frame(genelist)
qtl=as.data.frame(qtl_salt)
qtldist(geneset, genelist, qtl, plot=TRUE)
}
\keyword{gene}
\keyword{geneset}
\keyword{genelist}
\keyword{qtl}
