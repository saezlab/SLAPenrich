\name{LUAD_CaseStudy_clinicalInfos}
\alias{LUAD_CaseStudy_clinicalInfos}
\docType{data}
\title{
    Clinical informations for a cohort of 188 lung adenocarcinoma patients
}
\description{
A named binary matrix where column names are clinical factors, and the row names sample identifiers. A non-zero entry in position \eqn{i,j} of this matrix indicates the for the \eqn{i-}sample the in the \eqn{j-}factor is positive.
This matrix summarizes some clinical informations for a cohort of 188 lung adenocarcinoma patients of a public available dataset (see source). Particularly the smoking status of the patient-samples (former smoker, current smoker, never smoked, not available) and the 
bronchioalveolar carcinoma type (mucinous and not-mucinous).
This dataset is paired with \code{\link{LUAD_CaseStudy}}, which summarises the somatic variants found in the same cohort of patients.
}
\format{
 A named binary matrix with sample identifiers as row names and clinical factor identifiers as column names:
 i.e.
 format:
 num [1:188, 1:6] 1 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:188] "16530" "16594" "16596" "16600" ...
  ..$ : chr [1:6] "SS_NotAvailable" "SS_Former" "SS_CurrentSmoker" "SS_Never" ...
}
\details{
The sample identifiers on the rows match the column names of the integer matrix in \code{\link{LUAD_CaseStudy}}.
}
\source{
The dataset from which this matrix was derived has been studied in Ding et al, 2008. The clinical informations used to assemble this matrix are in Supplementary Table 15 of this publication (available at \url{http://genome.wustl.edu/pub/supplemental/tsp_nature_2008/})
}
\seealso{
\code{\link{LUAD_CaseStudy}}
}
\references{
Ding L, Getz G, Wheeler DA, Mardis ER, McLellan MD, Cibulskis K, et al. Somatic mutations affect key pathways in lung adenocarcinoma. Nature. 2008;455:1069-75.
}
\keyword{datasets}
