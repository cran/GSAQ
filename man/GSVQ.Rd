\name{GSVQ}
\alias{GSVQ}
\title{Gene Set Validation with QTL using Hyper-geometric test without gene sampling model}
\description{The function computes ths statisical significance value (p-value) for gene set validation using hypergeometric test.}
\usage{GSVQ(geneset, genelist, qtl)}
\arguments{
  \item{geneset}{geneset is a vector of characters representing the names of genes/ gene ids selected from the whole gene list by using a gene selection method.}
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names): where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs.}
  \item{qtl}{qtl is a Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs.}
}
\value{The function returns the statisical significance value (p-value) from Hyper-geometric test for validation of the selected gene set with qtl data.}

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
GSVQ(geneset, genelist, qtl)
}
\keyword{gene}
\keyword{geneset}
\keyword{genelist}
\keyword{qtl}