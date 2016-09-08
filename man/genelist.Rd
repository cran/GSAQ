\name{genelist}
\alias{genelist}
\docType{data}
\title{
A list of genes of rice  
}
\description{
This data is in form of a 200 by 3 dataframe with genes/gene ids as rownames. The first column represents the chromosomal location of the genes (chromosome number). The second coloumn represents start position of the genes in terms of basepairs (bps) and the third coloumn represents end position of genes in terms of basepairs (bps) in their respective chromosomes. 
}
\usage{data("genelist")}
\format{
  A data frame with 200 rows as genes and the columns represent the chromosomal locations, start positions and end positions of respective genes.
  \describe{
    \item{\code{Chr}}{chr represents the chromosomal location of the genes}
    \item{\code{Start}}{start represents the start position of the genes in their respective chromosomes}
    \item{\code{End}}{End represents the end position of the genes in their respective chromosomes}
  }
}
\details{
The data is created by taking 200 genes from the large number of genes from NCBI GEO database. The genomic location of the genes on the rice genome are obtained from MSU Rice Genome Annotation (Osa1).}
\source{
Gene Expression Omnibus: NCBI gene expression and hybridization array data repository.ncbi.nlm.nih.gov/geo/.
Ouyang S, Zhu W, Hamilton J, Lin H, Campbell M, et al. (2007) The TIGR Rice Genome Annotation Resource: improvements and new features. Nucleic Acids Research 35.
}
\examples{
data(genelist)
}
\keyword{genes}
\keyword{chromosome}
\keyword{start position}
\keyword{end position}
