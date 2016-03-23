\name{SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection}
\alias{SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection}
\docType{data}
\title{
Collection of pahtway gene-sets from Pathway-Commons post-processed for redundancy reduction
}
\description{
A list containing pathway gene-sets from multiple public resources, downloaded from Pathway-Commons and post-processed to reduce
their overlaps (see details).
}
\format{
   A list containing the following items:
  \describe{
    \item{\code{PATHWAY}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons name, or multiple Pathway-Commons name joined
    (separated by '//'), for the \eqn{i}-th pathway gene-set (or gene-set resulting from merging multiple pathways, see details);}
    \item{\code{SOURCE}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons description of the source of the                               \eqn{i}-th pathway, or sources of multiple merged pathways;}
    \item{\code{UNIPROTID}}{A list in which the \eqn{i}-th element is a string vector containing the uniprot identifiers of the genes                                 belonging to the \eqn{i}-th pathway;}
    \item{\code{HGNC_SYMBOL}}{A list in which the \eqn{i}-th element is a string vector containing the official HUGO symbols of the genes                                 belonging to the \eqn{i}-th pathway or multiple merged pathways}
    \item{\code{Ngenes}}{An integer vector in which the \eqn{i}-th element is the number of genes belonging to the \eqn{i}-th pathway;}
    \item{\code{backGround}}{A string vector containing the HUGO symbols of all the genes belonging to at least one pathway}
    \item{\code{miniSOURCE}}{A string vector in which the \eqn{i}-th entry contains the name of the source of the                                                      \eqn{i}-th pathway (panther, humancyc, pid or reactome);}
    \item{\code{includesTP53}}{A boolean vector whose \eqn{i}-th is \code{TRUE} if the \eqn{i}-th pathway contains TP53.}
  }
  Please note that the name of this list is \code{PATHCOM_HUMAN}.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection)
## maybe str(SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection) ; plot(SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection) ...
}
\keyword{datasets}
