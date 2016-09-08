\name{totqtlhit}
\alias{totqtlhit}
\title{Computation of total number of qtl-hits found in the whole gene space}
\description{It enable to Compute the total number qtl-hits found in the whole gene space or in the micro-array chip}
\usage{totqtlhit(genelist, qtl)}
\arguments{
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs.}
  \item{qtl}{qtl is a Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs.}
}
\value{The function returns a numeric value representing the total number of qtl-hits found in the whole gene list or in a micro-array chip.}

\author{Samarendra Das}
\examples{
data(genelist)
data(qtl_salt)
genelist=as.data.frame(genelist)
qtl=as.data.frame(qtl_salt)
totqtlhit(genelist, qtl)
}
\keyword{genelist}
\keyword{qtl}
\keyword{qtlhit}