\name{SLAPE.PATHCOM_HUMAN_nr_i_hu_2014}
\alias{SLAPE.PATHCOM_HUMAN_nr_i_hu_2014}
\docType{data}
\title{
Collection of pahtway gene-sets from Pathway-Commons (v2014) post-processed for redundancy reduction and to update composing gene name to recent HUGO gene symbols
}
\description{
A list containing pathway gene-sets from multiple public resources, downloaded from Pathway-Commons and post-processed to reduce
their overlaps (see details) and update gene names.
}
\format{
   A list containing the following items:
  \describe{
    \item{\code{PATHWAY}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons name, or multiple Pathway-Commons name joined
    (separated by '//'), for the \eqn{i}-th pathway gene-set (or gene-set resulting from merging multiple pathways, see details);}
    \item{\code{SOURCE}}{A string vector in which the \eqn{i}-th entry contains the Pathway-Commons description of the source of the                               \eqn{i}-th pathway, or sources of multiple merged pathways;}
    \item{\code{UNIPROTID}}{A list in which the \eqn{i}-th element is a string vector containing the uniprot identifiers of the genes                                 belonging to the \eqn{i}-th pathway;}
    \item{\code{HGNC_SYMBOL}}{A list in which the \eqn{i}-th element is a string vector containing the official HUGO symbols of the genes                                 belonging to the \eqn{i}-th pathway or multiple merged pathways, differently from
    \code{SLAPE.PATHCOM_HUMAN_nr_i_hu} in this object these symbols are updated to recent nomenclature;}
    \item{\code{Ngenes}}{An integer vector in which the \eqn{i}-th element is the number of genes belonging to the \eqn{i}-th pathway;}
    \item{\code{backGround}}{A string vector containing the HUGO symbols of all the genes belonging to at least one pathway;}
    \item{\code{miniSOURCE}}{A string vector in which the \eqn{i}-th entry contains the name of the source of the                                                      \eqn{i}-th pathway (panther, humancyc, pid or reactome);}
    \item{\code{includesTP53}}{A boolean vector whose \eqn{i}-th is \code{TRUE} if the \eqn{i}-th pathway contains TP53.}
  }
  Please note that the name of this list is \code{PATHCOM_HUMAN}.
}
\details{
This object was assembled from a collection of pathway gene sets from the Pathway Commons data portal. From this collection gene sets containing less than 4 genes were discarded. Additionally, in order to remove redundancies those gene sets i) corresponding to the same pathway across different resources or ii) with a large overlap (Jaccard index (\eqn{J}) \eqn{> 0.8}, as detailed below) were merged together by intersecting them. The gene sets resulting from these compressions were then added to the collection (with a joint pathway label) and those participating in at least one of these merging were discarded.
The final collection resulting from this pre-processing is composed by 1,636 gene sets, for a total amount of 8,056 unique genes included in at least one gene set.
Given two gene sets \eqn{P_1} and \eqn{P_2} the corresponding \eqn{J(P_1,P_2)} is defined as:

\deqn{J(P_1,P_2)=\frac{|P_1 \cap P_2|}{|P_1 \cup P_2|}}{J(P_1,P_2)=(|P_1 ∩ P_2|)/(|P_1 ∪ P_2|)}.

Additionally, all the pathway gene sets contained in this object are updated to recent official HUGO gene nomenclatures, using the informations
contained in the \code{SLAPE.hgnc.table} data object (which can be itself updated using the dedicated function \code{SLAPE.update_HGNC_Table}).

}
\source{
This list was assembled from the collection of pathway gene sets from the Pathway-Commons data portal (v4-201311) (Cerami et al, 2011) (\url{http://www.pathwaycommons.org/archives/PC2/v4/}). 
}
\references{
Cerami EG, Gross BE, Demir E, Rodchenkov I, Babur O, Anwar N, et al. Pathway Commons, a web resource for biological pathway data. Nucleic Acids Res. 2011;39:D685-90
}
\seealso{
\code{\link{SLAPE.PATHCOM_HUMAN}},
\code{\link{SLAPE.update_HGNC_Table}}
}
\keyword{datasets}
