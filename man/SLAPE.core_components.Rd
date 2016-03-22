\name{SLAPE.core_components}
\alias{SLAPE.core_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identification of core-component genes shared by multiple SLAPenriched pathways
}
\description{
This function identifies group of genes that are recurrently altered in the analysed dataset and that are shared by
multiple SLAPenriched pathways, thus are putatively leading the enrichment scores. Additionally this function generates pdf files
containing \emph{pathway-membership heatmaps} showing to which pathway each of the genes in the core-component belongs to,
together with barplots with alteration frequencies for all the genes in the core-components. Results are also stored in individual
Robjects.
}
\usage{
SLAPE.core_components(PFP, EM, PATH = "./", fdrth = Inf, exclcovth = 0,
                      PATH_COLLECTION)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PFP}{
    A list containing the SLAPenrich analysis results outputted by the \code{\link{SLAPE.analyse}} function while analysing
        the genomic dataset summarised by the \code{EM} genomic event matrix.
}
  \item{EM}{
    A sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the column names are sample identifiers, and     the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the presence of a somatic      mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene. This matrix must be the same that has been inputted to the
       \code{\link{SLAPE.analyse}} function to produce the results in the \code{PFP} list.
}
  \item{PATH}{
    String specifiying the path of the directory where the pdf file shoud be created.
}
  \item{fdrth}{
    The false discovery rate threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage).
      By default in the \code{PFP} with an FDR < 5\% are selected.
}
  \item{exclcovth}{
    The mutual exclusivity coverage threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage).     By default in the \code{PFP} with an exclusive coverage > 50\% are selected.
}
  \item{PATH_COLLECTION}{
The pathway collection that has been tested against the \code{EM} dataset with the \code{\link{SLAPE.analyse}} function to produce           the \code{PFP} list of results.
}
}
\details{
To identify shared core-components across significantly enriched pathways, the set of enriched pathways and their composing genes are modeled as a bipartite network, in which the first set of nodes contains one element per enriched pathway and the second set contains one element per each of the genes that are included in at least one of the enriched pathways.

In this network, a pathway node is connected with an edge to each of its composing gene nodes.

The resulting bipartite network is then mined for communities, i.e. groups of densely interconnected nodes by using a fast community detection algorithm based on a greedy strategy (Newman, 2004).

The resulting communities are finally saved into pdf files containing heatmaps where nodes in the first set (pathways) are on the columns by columns, nodes in the second set (genes) are on the rows and a not-empty cell in position i,j indicates that the i-th gene belongs to the j-th pathway.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Newman MEJ. Fast algorithm for detecting community structure in networks. Phys Rev E Stat Nonlin Soft Matter Phys. 2004;69:066133.

Csardi G, Nepusz T. The igraph software package for complex network research. InterJournal, Complex Systems. 2006;1695:38
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
This function makes use of the \code{fastgreedy.community} function of the \code{igraph} R package (Csardi and Nepusz, 2006).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.analyse}}
}
\examples{
#Loading the Genomic Event data object derived from variants annotations
#identified in 188 Lung Adenocarcinoma patients (Ding et al, 2008)
data(LUAD_CaseStudy_updatedGS)

#Loading KEGG pathway gene-set collection data object
data(SLAPE.20160211_MSigDB_KEGG_hugoUpdated)

#Loading genome-wide total exonic block lengths
data(SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209)

#Running SLAPenrich analysis
PFPw<-SLAPE.analyse(EM = LUAD_CaseStudy_ugs,
                    PATH_COLLECTION = KEGG_PATH,
                    BACKGROUNDpopulation = rownames(LUAD_CaseStudy_ugs),
                    GeneLenghts = GECOBLenghts)

#Generating pdf files containing heatmaps of the core-components
#of SLAPenriched pathway with an FDR < 5\% and exclusive coverage > 80\%.
#The pdf files are saved in the current working directory.
SLAPE.core_components(EM = LUAD_CaseStudy_ugs,PFP = PFPw,
                      exclcovth = 80,fdrth = 5,
                      PATH = './',PATH_COLLECTION = KEGG_PATH)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Analysis}
\keyword{Result-exploration}