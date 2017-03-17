\name{LUAD_CaseStudy}
\alias{LUAD_CaseStudy}
\docType{data}
\title{
    Genomic event matrix derived from variants found in a cohort of 188 lung adenocarcinoma patients
}
\description{
A sparse integer matrix where column names are sample identifiers, and the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the number of somatic point mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene.
This matrix summarizes the somatic variants of a cohort of 188 lung adenocarcinoma patients of a public available dataset (see source).
}
\format{
  A named integer matrix with HUGO official gene symbols as row names and sample identifiers as column names:
  i.e.
  format: num [1:356, 1:163] 1 0 0 0 0 0 0 0 0 0 ...
 
 - attr(*, "dimnames")=List of 2
 
  ..$ : chr [1:356] "ABL1" "ABL2" "ACVR1B" "ACVR1C" ...
  ..$ : chr [1:163] "16770" "16646" "17741" "16915" ...

}
\source{
The dataset from which this matrix was derived has been studied in Ding et al, 2008. The variant annotations used to assemble this matrix are in Supplementary Table 2 of this publication (available at \url{http://genome.wustl.edu/pub/supplemental/tsp_nature_2008/})
}
\references{
Ding L, Getz G, Wheeler DA, Mardis ER, McLellan MD, Cibulskis K, et al. Somatic mutations affect key pathways in lung adenocarcinoma. Nature. 2008;455:1069-75.
}
\seealso{
\code{\link{LUAD_CaseStudy_ugs}}
}
\keyword{datasets}
