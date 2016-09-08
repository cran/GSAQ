\name{qtl_salt}
\alias{qtl_salt}
\docType{data}
\title{
A list of salt responsive Quantitative Trait Loci of rice  
}
\description{
This data is in form of a 13 by 3 dataframe with qtls/qtl ids as rownames. The first column reoresents the chromosomal location of the respective qtls (chromosome number). The second coloumn represents start position of the qtls in terms of basepairs (bps) and the third coloumn represents end position of qtls in terms of basepairs (bps) in their respective chromosomes. 
}
\usage{data("qtl_salt")}
\format{
  A data frame with 13 rows as qtl and the columns represent the chromosomal locations, start positions and end positions of respective qtls.
  \describe{
    \item{\code{Chr}}{chr represents the chromosomal location of the qtls}
    \item{\code{Start}}{start represents the start position of the qtls in their respective chromosomes}
    \item{\code{End}}{End represents the end position of the qtls in their respective chromosomes}
  }
}
\details{
The data is created by taking 13 unique salt responsive qtls from the Gramene QTL database. The genomic locations of these QTLs on rice genome are obtained using Gramene annotation of MSU Rice Genome Annotation (Osa1).}
\source{
Gramene QTL library (http://www.gramene.org/qtl/).
Ouyang S, Zhu W, Hamilton J, Lin H, Campbell M, et al. (2007) The TIGR Rice Genome Annotation Resource: improvements and new features. Nucleic Acids Research 35: D883-D887.
}
\examples{
data(qtl_salt)
}
\keyword{qtls}
\keyword{chromosome}
\keyword{start position}
\keyword{end position}
\keyword{base pair}
