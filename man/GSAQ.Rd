\name{GSAQ}
\alias{GSAQ}
\title{Gene Set Analysis with Quantitative Trait Loci with gene sampling model}
\description{The function computes the statistical significance value (p-value) from gene set analysis test with QTL for the test H0: Genes in the selected geneset are at most as often overlapped with the QTL regions as the genes in not selected geneset; against H1: Genes in the geneset are more often overlapped with the QTL regions as compared to genes in not selected geneset.}
\usage{GSAQ(geneset, genelist, qtl, SampleSize, K, method)}
\arguments{
  \item{geneset}{geneset is a vector of characters representing the names of genes/ gene ids selected from the whole gene list by using a gene selection method.}
  \item{genelist}{genelist is a N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene space: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.}
  \item{qtl}{qtl is a Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.}
  \item{SampleSize}{SampleSize is a numeric constant representing the size of the gene sample drawn from the geneset using the gene sampling model (SampleSize must be less than the size of geneset).}
  \item{K}{K is a numeric constant representing the number of gene samples of size equal to SampleSize will be drawn by the using gene sampling model.}
  \item{method}{method is a character string indicating which method for final p-value (combining p-values for various gene samples) is to be computed. One of "meanp", "sump", "logit", "sumz"or "logp" (default) can be abbreviated and used.}
}
\value{The function returns the final statistical significance value (p-value) from Gene set Analysis with QTL test.}
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
GSAQ(geneset, genelist, qtl, SampleSize=30, K=50, method="meanp")
}
\keyword{gene}
\keyword{geneset}
\keyword{genelist}
\keyword{qtl}
\keyword{p value}
\keyword{gene sampling model}