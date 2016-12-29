\name{SLAPE.serialPathVis}
\alias{SLAPE.serialPathVis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Systematic generation of heatmaps of the alteration matrices for SLAPenriched pathways and barplots with statistical scores
}
\description{
    This function generates pdf files containing heatmaps of the alteration matrices for SLAPenriched pathways (with user-defined enrichment FDR and exclusive coverage) across the samples of the analysed dataset, after permuting rows and columns to highlight trend of mutual exclusivity in the alteration-patterns.
    
    Additionally, it generates, a pdf files with three barplots indicating, for each pathway, respectively: (i) the number of mutated genes across samples; (ii) the probabilities of the pathway under consideration to be altered across samples, together with the expected number of samples with alteration in the pathway under consideration; (iii) The observed pathway alteration status across samples, together with the observed total number of samples with alteration in the pathway under consideration.
}
\usage{
SLAPE.serialPathVis(EM, PFP, fdrth = 5, exCovTh = 50, PATH = "./", PATH_COLLECTION)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EM}{
A sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the column names are sample identifiers, and     the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the presence of a somatic      mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene. This matrix must be the same that has been inputted to the
       \code{\link{SLAPE.analyse}} function to produce the results in the \code{PFP} list.
}
  \item{PFP}{
A list containing the SLAPenrich analysis results outputted by the \code{\link{SLAPE.analyse}} function while analysing
        the genomic dataset summarised by the \code{EM} genomic event matrix.
}
  \item{fdrth}{
The false discovery rate threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage).
      By default in the \code{PFP} with an FDR < 5\% are selected.
}
  \item{exCovTh}{
    The mutual exclusivity coverage threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage).     By default in the \code{PFP} with an exclusive coverage > 50\% are selected.
}
  \item{PATH}{
    String specifiying the path of the directory where the pdf file shoud be created.
}
  \item{PATH_COLLECTION}{
The pathway collection that has been tested against the \code{EM} dataset with the \code{\link{SLAPE.analyse}} function to produce           the \code{PFP} list of results.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
This function calls iteratively the \code{\link{SLAPE.pathvis}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.analyse}}, \code{\link{SLAPE.pathvis}}
}
\examples{
    #Loading the Genomic Event data object derived from variants annotations
    #identified in 188 Lung Adenocarcinoma patients (Ding et al, 2008)
    data(LUAD_CaseStudy_updatedGS)

    #Loading KEGG pathway gene-set collection data object
    data(SLAPE.MSigDB_KEGG_hugoUpdated)

    #Loading genome-wide total exonic block lengths
    data(SLAPE.all_genes_exonic_content_block_lengths_ensemble)

    #Running SLAPenrich analysis
    PFPw<-SLAPE.analyse(EM = LUAD_CaseStudy_ugs,
                    PATH_COLLECTION = KEGG_PATH,
                    BACKGROUNDpopulation = rownames(LUAD_CaseStudy_ugs),
                    GeneLenghts = GECOBLenghts)

    #Generating pdf files containing heatmaps of the mutual-exclusivity
    #sorted pathway alteration matrices, for an SLAPenriched pathway with
    #an exclusive coverage > 80\%, and barplots with statistical scores.
    #The pdf files are saved in the current working directory.
    SLAPE.serialPathVis(EM = LUAD_CaseStudy_ugs,PFP = PFPw,
                    exCovTh = 80,fdrth = 5,
                    PATH = './',PATH_COLLECTION = KEGG_PATH)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Result-exploration}
