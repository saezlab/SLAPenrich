\name{SLAPE.20160211_MSigDB_KEGG_hugoUpdated}
\alias{SLAPE.20160211_MSigDB_KEGG_hugoUpdated}
\docType{data}
\title{
Collection of KEGG pahtway gene-sets from the Molecular Signature Database post-processed to update composing gene name to recent HUGO
gene symbols
}
\description{
A list containing KEGG pathway gene-sets downloaded from the Molecular Signature Database (Subramanian et al, 2005) and post-processed to update gene names to recent HUGO nomenclature.
}
\format{
   A list containing the following items:
  \describe{
    \item{\code{PATHWAY}}{A string vector in which the \eqn{i}-th entry contains the KEGG pathway name for the \eqn{i}-th pathway gene-set;}
    \item{\code{SOURCE}}{A string vector in which the \eqn{i}-th entry contains the source of the \eqn{i}-th pathway;}
    \item{\code{HGNC_SYMBOL}}{A list in which the \eqn{i}-th element is a string vector containing the official HUGO symbols of the genes                                 belonging to the \eqn{i}-th pathway or multiple merged pathways, updated to recent nomenclature;}
    \item{\code{Ngenes}}{An integer vector in which the \eqn{i}-th element is the number of genes belonging to the \eqn{i}-th pathway;}
    \item{\code{backGround}}{A string vector containing the HUGO symbols of all the genes belonging to at least one pathway;}
}
  Please note that the name of this list is \code{KEGG_PATH}.
}
\details{
All the pathway gene sets contained in this object are updated to recent official HUGO gene nomenclatures, using the informations
contained in the \code{SLAPE.hgnc.table_20160210} data object (which can be itself updated using the dedicated function \code{SLAPE.update_HGNC_Table}
}
\source{
This list was assembled from a collection of KEGG pathway gene sets from the Molecular Signature Databse (Subramanian et al, 2005) (\url{http://http://software.broadinstitute.org/gsea/msigdb}). 
}
\references{
Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP.
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
Proc Natl Acad Sci U S A. 2005 Oct 25;102(43):15545-50
}
\keyword{datasets}
