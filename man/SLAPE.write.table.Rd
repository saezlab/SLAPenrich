\name{SLAPE.write.table}
\alias{SLAPE.write.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Writing SLAPenrich results in a csv file
}
\description{
This function takes in input a list of results outputted by the \code{\link{SLAPE.analyse}}, selects enriched pathways
according to user-defined significance and mutual exclusivity coverage thresholds and creates an easy to explore csv file
with selected enriched pathways.
}
\usage{
SLAPE.write.table(PFP,
                  EM,
                  filename = "",
                  fdrth = Inf,
                  exclcovth = 0,
                  PATH_COLLECTION,
                  GeneLenghts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PFP}{
    A list containing the SLAPenrich analysis results outputted by the \code{\link{SLAPE.analyse}} function while analysing
    the genomic dataset summarised by the \code{EM} genomic event matrix.
}
  \item{EM}{
       A sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the column names are sample identifiers, and the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the presence of a somatic mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene. This matrix must be the same that has been inputted to the
       \code{\link{SLAPE.analyse}} function to produce the results in the \code{PFP} list.
}
  \item{filename}{
  String specifiying the path of the csv file to be created (including its name).
}
  \item{fdrth}{
      The false discovery rate threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage).
      By default all the pathway included in the \code{PFP} list are selected.
}
  \item{exclcovth}{
      The mutual exclusivity coverage threshold that should be used to select SLAPenriched pathways from the \code{PFP} list (percentage). By default all the pathway included in the \code{PFP} list are selected.
}
  \item{PATH_COLLECTION}{
        The pathway collection that has been tested against the \code{EM} dataset with the \code{\link{SLAPE.analyse}} function to produce the \code{PFP} list of results.
}
  \item{GeneLenghts}{
    The named vector containing the genome-wide total exonic block lengths that has been used by the \code{\link{SLAPE.analyse}} function to produce the \code{PFP} list of results.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
    \code{\link{SLAPE.analyse}}
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
                  
    #Selecting pathway enriched at a 5\% FDR,
    #that have a 50\% mutual exclusivity coverage and writing them
    #in a csv file
    SLAPE.write.table(PFP = PFPw,
                      EM = LUAD_CaseStudy_ugs,
                      filename = "./LungDS_KEGG_enrichments.csv",
                      fdrth=5, exclcovth = 50, PATH_COLLECTION = KEGG_PATH,
                      GeneLenghts = GECOBLenghts)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Result-exploration}