\name{geneqtl}
\alias{geneqtl}
\title{List of the selected genes along with their corresponding overlapped QTL}
\description{The function enables to obtain list of the selected genes along with the corresponding overlapped Quantitative Trait Loci (QTL) ids/names along with their genomic positions.}
\usage{geneqtl(geneset, genelist, qtl)}
\arguments{
  \item{geneset}{geneset is a vector of characters representing the names of genes/ gene ids selected from the whole gene list/space by using a gene selection method.}
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.}
  \item{qtl}{qtl is a Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.}
}
\value{The function returns a list with two components. First component returns the list of selected genes along with their overlapped QTL ids/names. Second component gives the list of selected genes with their overlapped QTL ids/names and their respective genomic positions.}
\author{Samarendra Das}
\examples{
data(rice_salt)
data(genelist)
data(qtl_salt)
x=as.data.frame(rice_salt[-1,])
y=as.numeric(rice_salt[1,])
genelist=as.data.frame(genelist)
qtl=as.data.frame(qtl_salt)
geneset= GeneSelect(x, y, s=50, method="t-score")$selectgenes
geneqtl(geneset, genelist, qtl)

}
\keyword{gene}
\keyword{geneset}
\keyword{genelist}
\keyword{QTL}