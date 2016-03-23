\name{SLAPE.diff_SLAPE_analysis}
\alias{SLAPE.diff_SLAPE_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential SLAPenrichment analysis
}
\description{
This function allows the identification of pathways that are differentially enriched (thus frequently altered) across two sub-populations of samples of the same input dataset. Similarly to differential gene expression analysis, the two sub-populations to be contrasted are defined through a contrast matrix. 
}
\usage{
SLAPE.diff_SLAPE_analysis(EM, contrastMatrix,
                          positiveCondition, negativeCondition,
                          show_progress = TRUE, display = TRUE,
                          correctionMethod = "fdr",
                          path_probability = "Bernoulli",
                          NSAMPLES = 1, NGENES = 1,
                          accExLength = TRUE,
                          BACKGROUNDpopulation = NULL,
                          GeneLenghts,
                          PATH_COLLECTION,
                          SLAPE.FDRth = 5, PATH = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EM}{A sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the column names are sample identifiers, and the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the presence of a somatic mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene. If the matrix contains integer entries then these values are deemed as the number of somatic point mutations harbored by a given sample in a given gene (these values will be considered if the analysis takes into account of the gene exonic lengths, or converted in binary values otherwise).
}
  \item{contrastMatrix}{
        A binary matrix specifying which sample is included in which sub-population. The row names of this matrix are sample identifiers (and must match the column names of the \code{EM} dataset). The column names indicate the sub-population identifiers. A 1 in the position i,j of such a matrix indicates that the i-th sample is included in the sub-population corresponding to the j-th condition.
}
  \item{positiveCondition}{
    String indicating one of the two sub-populations of samples to be contrasted (the positive population). It should match a column header of the           \code{contrastMatrix}.
}
  \item{negativeCondition}{
  String indicating one of the two sub-populations of samples be contrasted (the negative population). It should match a column header of the           \code{contrastMatrix}.
}
  \item{show_progress}{
        Boolean parameter determining if a progress bar should be visualized during the execution of the analysis (default) or not.
}
  \item{display}{
        Boolean parameter determining if result figures should be displayed and saved.
}
  \item{correctionMethod}{
    A string indicating which method should be used to correct pathway enrichment p-values for multiple hypothesis testing, in the
    two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls). Possible values for this parameter are all the values for the     \code{method} parameter of the R \code{p.adjust} function, plus ``qvalue'' for the Storey -Tibshirani method (Storey and Tibshirani, 2003).
}
  \item{path_probability}{
    A string specifying which model should be used to compute the sample-wise pathway alteration probabilities, in the
    two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls). Possible values for this parameter are ``Bernoulli'' (default)     and ''HyperGeom''.
}
  \item{NSAMPLES}{
    The minimal number of samples of \code{EM} in which at least one gene belonging to a given pathway should be mutated in order for that        pathway to be tested for alteration enrichments at the sample population level, in the
    two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls).
}
  \item{NGENES}{
    The minimal number of genes of a given pathway \eqn{P} that must be mutated in at least one sample of the \code{EM} in order for that         pathway to be tested for alteration enrichments at the sample population level, in the
    two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls).
}
  \item{accExLength}{
    Boolean parameter determining whether the sample-wise pathway alteration probability model in the
    two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls) should take into account of the total exonic block lengths of the genes in the pathways and analyzed dataset (default) or not (see details), when using an Hypergeometric model(as specified by the \code{path_probability} parameter). This parameter is neglected if a Bernoulli model is used for these probabilities instead of a Hypergeometric model.
}
  \item{BACKGROUNDpopulation}{
 A string vector containing the official HUGO symbols of the gene background population used to compute the sample-wise pathway alteration probabilities in the two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls). If \code{NULL} (default) then the population of all the genes included in at least one pathway of the collection specified in the \code{PATH_COLLECTION} parameter
}
  \item{GeneLenghts}{
    A named vector containing the genome-wide total exonic block lengths to be used in the two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls).Names of this vector are official HUGO gene symbol. This is available in the
        \code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209 object}}.
        An updated version of this vector can be assembled using the
        \code{\link{SLAPE.compute_gene_exon_content_block_lengths}} function.
}
  \item{PATH_COLLECTION}{
    A list containing the pathway collection to be tested on the \code{EM} dataset for alteration enrichments at the sample population level,
     in the two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls). Several collections are included in the package as data      object.
        See for example
        \code{\link{SLAPE.20140608_PATHCOM_HUMAN}}
        or
        \code{\link{SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated}}
        for a description of the fields required in this list.
}
  \item{SLAPE.FDRth}{
        The false discovery rate threshold to be considered for selecting significant pathways in at least one of the two individual SLAPenrich analyses (therefore \code{SLAPE.analyse} calls), and for which differential enrichment scores should be computed.
}
  \item{PATH}{
        String specifiying the path of the directory where the pdf file containing results shoud be created.
}
}
\details{
This function first performs two individual SLAPenrichment analyses using the \code{SLAPE.analyse} function, on the user-defined two sub-populations of samples, yielding two lists of results. The pathways that are significantly enriched in at least one of the two result lists (according to a user defined false discovery rate (FDR) threshold) are then selected and, for each of them, a differential enrichment score is computed as:
\deqn{\Delta_{A,B}(P) = -\log_{10}\mathrm{FDR}_A(P)+\log_{10}\mathrm{FDR}_B(P),}{∆_(A,B) (P)= -log_10 FDR_A(P)+log_10 FDR_B(P),}
where \eqn{A} and \eqn{B} are the two contrasted sub-populations (respectively, positive and negative) and \eqn{\mathrm{FDR}_A}{FDR_A} and \eqn{\mathrm{FDR}_B}{FDR_B} are the two SLAPenrichment FDRs obtained in the two corresponding individual analyses. 

Results are visualised at the level of the inputted alterations across the two contrasted population, on the domain of the differentially enriched pathways as well has heatmaps and barplots of the differential enrichment scores.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Storey JD, Tibshirani R. Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences of the United States of America. 2003;100:9440–5.
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Loading the Genomic Event data object derived from variants annotations
#identified in 188 Lung Adenocarcinoma patients (Ding et al, 2008)
data(LUAD_CaseStudy_updatedGS)

#Loading KEGG pathway gene-set collection data object
data(
    SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated)

#Loading genome-wide total exonic block lengths
data(SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209)

#Loading clinical infos for 188 Lung Adenocarcinoma patients
#(Ding et al, 2008)
data(LUAD_CaseStudy_clinicalInfos)

#Performing differential SLAPenrichment analysis comparing
#Smokers Vs. Non Smokers. Pdf files with result figures are saved in
#the current working directory
RES<-
    SLAPE.diff_SLAPE_analysis(EM = LUAD_CaseStudy_ugs,contrastMatrix = LUAD_CaseStudy_clinicalInfos,
                              BACKGROUNDpopulation = rownames(LUAD_CaseStudy_ugs),
                              SLAPE.FDRth = 5,display = TRUE,
                              positiveCondition = 'SS_CurrentSmoker',
                              negativeCondition = 'SS_Never',
                              PATH_COLLECTION = PATHCOM_HUMAN,
                              GeneLenghts = GECOBLenghts,
                              PATH = './')

#Showing the top-10 most differentially enriched pathways in the Smokers population
RES[1:10,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Analysis}
\keyword{Result-exploration}
