\name{SLAPE.PATHCOM_HUMAN}
\alias{SLAPE.PATHCOM_HUMAN}
\docType{data}
\title{
Collection of pahtway gene-sets from Pathway-Commons
}
\description{
A list containing pathway gene-sets from multiple public resources, downloaded from Pathway-Commons.
}
\format{
  A list containing the following items:
  \describe{
    \item{\code{PATHWAY}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons name of the \eqn{i}-th pathway;}
    \item{\code{SOURCE}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons description of the source of the                               \eqn{i}-th pathway;}
    \item{\code{UNIPROTID}}{A list in which the \eqn{i}-th element is a string vector containing the uniprot identifiers of the genes                                 belonging to the \eqn{i}-th pathway;}
    \item{\code{HGNC_SYMBOL}}{A list in which the \eqn{i}-th element is a string vector containing the official HUGO symbols of the genes                                 belonging to the \eqn{i}-th pathway;}
    \item{\code{Ngenes}}{An integer vector in which the \eqn{i}-th element is the number of genes belonging to the \eqn{i}-th pathway;}
    \item{\code{backGround}}{A string vector containing the HUGO symbols of all the genes belonging to at least one pathway}
    \item{\code{miniSOURCE}}{A string vector in which the \eqn{i}-th entry contains the name of the source of the                                                      \eqn{i}-th pathway (panther, humancyc, pid or reactome);}
    \item{\code{includesTP53}}{A boolean vector whose \eqn{i}-th is \code{TRUE} if the \eqn{i}-th pathway contains TP53.}
  }
  Please note that the name of this list is \code{PATHCOM_HUMAN}.
}
\source{
    This list was assembled from the collection of pathway gene sets from the Pathway-Commons data portal (v4-201311) (Cerami et al, 2011) (\url{http://www.pathwaycommons.org/archives/PC2/v4-201311/}). 
}
\references{
Cerami EG, Gross BE, Demir E, Rodchenkov I, Babur O, Anwar N, et al. Pathway Commons, a web resource for biological pathway data. Nucleic Acids Res. 2011;39:D685-90
}

\keyword{datasets}
