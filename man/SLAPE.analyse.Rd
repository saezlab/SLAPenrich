
\name{SLAPE.analyse}
\alias{SLAPE.analyse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing Sample population Level Analysis of Pathway alteration Enrichments (SLAPenrich).
}
\description{
}
\usage{
SLAPE.analyse(EM,
              show_progress = TRUE,
              correctionMethod = "fdr",
              NSAMPLES = 1,
              NGENES = 1,
              accExLength = TRUE,
              BACKGROUNDpopulation = NULL,
              PATH_COLLECTION,
              path_probability = "Bernoulli",
              GeneLenghts,
              RHO=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EM}{A sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the column names are sample identifiers, and the row names official HUGO gene symbols. A non-zero entry in position \eqn{i,j} of this matrix indicates the presence of a somatic mutations harbored by the \eqn{j-}sample in the \eqn{i-}gene. If the matrix contains integer entries then these values are deemed as the number of somatic point mutations harbored by a given sample in a given gene (these values will be considered if the analysis takes into account of the gene exonic lengths, or converted in binary values otherwise).
}
  \item{show_progress}{
    Boolean parameter determining if a progress bar should be visualized during the execution of the analysis (default) or not.
}
  \item{correctionMethod}{
    A string indicating which method should be used to correct pathway enrichment p-values for multiple hypothesis testing. Possible values for this parameter are all the values for the \code{method} parameter of the R \code{p.adjust} function, plus ``qvalue'' for the Storey-Tibshirani method (Storey and Tibshirani, 2003).
}
  \item{NSAMPLES}{
    The minimal number of samples of \code{EM} in which at least one gene belonging to a given pathway should be mutated in order for that pathway to be tested for alteration enrichments at the sample population level.
}
  \item{NGENES}{
    The minimal number of genes of a given pathway \eqn{P} that must be mutated in at least one sample of the \code{EM} in order for that pathway to be tested for alteration enrichments at the sample population level.
}
  \item{accExLength}{
        Boolean parameter determining whether the sample-wise pathway alteration probability model should take into account of the total exonic block lengths of the genes in the pathways and analyzed dataset (default) or not (see details), when using an Hypergeometric model(as specified by the \code{path_probability} parameter). This parameter is neglected if a Bernoulli model is used for these probabilities instead of a Hypergeometric model.
}
  \item{BACKGROUNDpopulation}{
        A string vector containing the official HUGO symbols of the gene background population used to compute the sample-wise pathway alteration probabilities (see details). If \code{NULL} (default) then the population of all the genes included in at least one pathway of the collection specified in the \code{PATH_COLLECTION} parameter
}
  \item{PATH_COLLECTION}{
        A list containing the pathway collection to be tested on the \code{EM} dataset for alteration enrichments at the sample population level. Several collections are included in the package as data object.
        See for example
        \code{\link{SLAPE.PATHCOM_HUMAN}}
        or
        \code{\link{SLAPE.PATHCOM_HUMAN_nr_i_hu_2016}}
        for a description of the fields required in this list.
}
  \item{path_probability}{
        A string specifying which model should be used to compute the sample-wise pathway alteration probabilities (see details). Possible values for this parameter are ``Bernoulli'' (default) and ''HyperGeom''.
}
  \item{GeneLenghts}{
        A named vector containing the genome-wide total exonic block lengths. Names of this vector are official HUGO gene symbol. This is available in the
        \code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble}} object.
        An updated version of this vector can be assembled using the
        \code{\link{SLAPE.compute_gene_exon_content_block_lengths}} function.
}
  \item{RHO}{
        The mutation rate to be used to compute the sample-wise pathway alteration probabilities with the Bernoulli model. If \code{NA} (default) then this value is estimated directly by the \code{EM} dataset. The value of this parameter is ignored if an Hypergeometric model is used to compute these probabilities (this is specified by the \code{path_probability} parameter).
}
}
\details{
This is the core analysis function implementing the statistical framework described in (Iorio et al, in preparation).

It takes in input a dataset (the parameter \code{EM}) stored in a sparse binary matrix, or a sparse matrix with integer non-null entries. In this matrix the columns correspond to samples, the rows correspond to genes and a non-zero entry indicate the presence of a somatic mutations harbored by a given sample in a given gene. If the matrix contains integer entries then they are deemed as the number of somatic point mutations harbored by a given sample in a given gene (these values will be considered if the analysis takes into account of the gene exonic lengths, or converted in binary values otherwise, see below).

For each pathway gene-set \eqn{P} in the pathway collection specified in the parameter \code{PATH_COLLECTION}, and an inputted genomic dataset (in \code{EM}), this function computes first of all a vector of probabilities \eqn{\pi = \{p_i\}}{π = {p_i}}
quantifying how likely each sample is to harbor at least one somatic mutation in a gene belonging to \eqn{P}, by random chance.

These probabilities are computed (as specified by the parameter \code{path_probability}) by default using a Bernoulli model accounting for the total exonic block lengths of all the genes belonging to \eqn{P}, and the expected or observed background mutation rate \eqn{\rho}{ρ}:

\deqn{p_i = \mathrm{Pr}(X_i \geq 1) = 1 - \mathrm{exp} (-\rho k' )}{p_i= Pr(X_i ≥ 1) = 1 - exp (-ρk' )}

where \eqn{k' = \sum_{g \in P}{l(g)}}{k' = ∑_(g∈P) l(g)} is the sum of the exonic block length of all the genes \eqn{g} in pathway \eqn{P}, and \eqn{\rho}{ρ} the background mutation rate, which can be estimated from the input dataset directly or set to established estimated values (such as \eqn{10^{-6}/\mathrm{nucleotide}}{10^-6/nucleotide}) (Greenman et al, 2007; Stratton et al, 2009; Ding et al, 2008, Sjoblom et al, 2006), depending on the parameter \code{RHO} being set to its default \code{NA} value or not, and \eqn{X_i} is the total number of genes belonging to \eqn{P} that are mutated in sample \eqn{s_i}.

Alternatively, these probabilities can be computed through a complementary cumulative hypergeometric distribution evaluated at \eqn{X = 0}, and taking into account of the mutation burden of the samples \eqn{n_i}, the size of \eqn{P} in terms of number of genes \eqn{k}
the number of genes in the background-population \eqn{N}
\deqn{p_i = \mathrm{Pr}(X_i \geq 1) = \sum_{x=1}^{k} H(x,N,k,n_i)}{p_i = Pr(X_i ≥ 1) = ∑_(x=1)^k H(x,N,k,n_i)}

where \eqn{H} is the probability mass function of a hypergeometric distribution with parameters \eqn{x}, \eqn{N}, \eqn{k} and \eqn{n_i}.
In this case the parameter \code{accExLength} must be set to \code{FALSE}.

To take into account of the total exonic block lengths of the genes in \eqn{P} and the background population both then the parameter \code{accExLength} must be set to its default parameter (\code{TRUE}). In this case the parameter of \eqn{H} will be:

\enumerate{
    \item \eqn{N' = \sum_{g \in G}{l(g)}}{N'=∑_(g∈G) l(g)} : the sum of all the exonic content block lengths of all the genes in the background population \eqn{G};
    \item \eqn{k' = \sum_{g \in P}{l(g)}}{k' = ∑_(g∈P) l(g)} : the sum of the exonic block lengths of all the genes in pathway \eqn{P})
    \item \eqn{n_i} : the total number of individual alterations involving genes of the pathway \eqn{P} in sample \eqn{s_i}).
  }
  
In each case the gene background population \eqn{G} can be defined by the user (through the parameter \code{BACKGROUNDpopulation}) or assembled pooling together all the genes belonging to at least one pathway of the collection specified in \code{PATH_COLLECTION}.

After the vector of probabilities \eqn{\pi}{π} has been computed, this function computes a pathway alteration score, quantifying the deviance of the number of samples in the datasets harbouring at least a somatic mutation in at least one gene of \eqn{P}, \eqn{O(P)}, from its random expectation \eqn{E(P)}:

\deqn{A(P) = \log_{10}{\frac{O(P)}{E(P)}},}{A(P) = log10{O(P)/E(P)},}

where \eqn{E(P)=\sum_{i}^m{p_i}}{E(P)= ∑_(i=1)^m p_i} and \eqn{m} is the number of samples in the anlyzed dataset.

Finally, this function computes an alteration score significance with a p-value against the null hypothesis: ``O(P) is drawn from a Poisson binomial distribution with \eqn{\pi = \{p_i\}}{π = {p_i}} success probabilities''. This comes from the observation that if there is no tendency for a given pathway to be recurrently mutated across \eqn{m} samples of the datasets then each of these samples can be considered as the observation of a single Bernoulli trial (in a series of \eqn{m} of them), where the event under consideration in the i-th trial is ``At least one gene belonging to \eqn{P} is mutated in the i-th sample''. The success probability of this event is given by \eqn{p_i}.

After alteration scores and corresponding significance have been assessed for all the pathways considered in the analysis, resulting p-values are corrected for multiple hypothesis testing with a user-defined method, specified by the parameter \code{correctionMethod}.

For each of the tested pathways an exclusive coverage score quantifying the tendency of the genes in a given pathway is also computed.
This is quantified as the ration between the number of samples where a unique gene belonging to a given pathway is altered divided by the number of samples where at least one gene belonging to given pathway is altered.

}
\value{

    A list containing the results of the sample-population level analysis of pathway alterations, containing the following items:

    \item{pathway_id}{A numerical vector containing one numerica value for each tested pathway. This is the index of the pathway in the vectors and lists of the pathway collection object specified in the \code{PATH_COLLECTION} input parameter;}
    
    \item{pathway_EM}{The pathway-level alteration matrix: a binary matrix with pathway indexes on the rows and sample identifiers on the columns, and binary entries indicating whether a given pathway is genomically altered in a given sample. The column names are the same of the \code{EM} input dataset;}


    \item{pathway_Probability}{The sample-wise pathway alteration probabilities. A numerical matrix with pathway indexes on the rows, sample identifiers on the columns where the generic entry \eqn{i,j} quantifies the likelyhood of the \eqn{i}-pathway to be genomically altered in the \eqn{j}-sample according to the model selected through the input parameter specifications (see Details and Arguments);}
    
    \item{pathway_mus}{The expected number of samples in which each pathway is altered: a named vector in which each entry specifies for each pathway (whose index is the entry name) the expected number of samples where it is altered according to the model selected through the input parameter specifications;}

    \item{pathway_logOddRatios}{The pathway alteration scores: a named vector in which each entry specifies for each pathway (whose index is the entry name) the \eqn{\log_{10}}{log_10} ratio between the number of samples in which that pathway is altered divided by the expected number of samples in which that pathway is altered (this is computed according to the model selected through the input parameter specifications);}

    \item{pathway_pvals}{Significance of the pathway alteration scores: a named vector of p-values against the null hypothesis that there is no association between the disease population summarized in the \code{EM} input dataset and a pathway P. The entry-names of this vector are the pathway indexes;}


    \item{pathway_perc_fdr}{Significance of the pathway alteration scores after correction for multiple hypotheses testing: a named vector of adjusted p-values (in the form of false discovery rate percentages) against the null hypothesis that there is no association between the disease population summarized in the \code{EM} input dataset and a pathway P. The entry-names of this vector are the pathway indexes.}

    \item{pathwayExclusive_coverage}{The mutual exclusive coverage of each pathway: a named vector containing the mutual exclusivity coverage score of each pathway (see details). The names of this vector are the pathway indexes.}
    
    \item{pathway_individualEMs}{The individual pathway alteration matrices. A named list with an entry for each tested pathway. Each of this entry is a binary matrix with all the genes belonging to the pathway under consideration on the rows, sample identifiers on the columns and binary entries specifing whether a given gene is altered in a given sample. The column names of these matrix are the same of the input \code{EM} dataset, whereas the entry-names of this list are the pathway indexes.}

}
\references{
Ding L, Getz G, Wheeler DA, Mardis ER, McLellan MD, Cibulskis K, et al. Somatic mutations affect key pathways in lung adenocarcinoma. Nature. 2008;455:1069-75.

Greenman C, Stephens P, Smith R, Dalgliesh GL, Hunter C, Bignell G, et al. Patterns of somatic mutation in human cancer genomes. Nature. 2007;446:153-8.

Iorio F, Garcia-Alonso L, Buendia JE, Brammeld J, Wille DR, McDermott U, Saez-Rodriguez J. Identification and visualization of population-level enrichments of pathway alterations in cancer genomes with SLAPenrich. [In preparation]


Sjoblom T, Jones S, Wood LD, Parsons DW, Lin J, Barber TD, et al. The consensus coding sequences of human breast and colorectal cancers. Science. 2006;314:268-74.

Stratton MR, Campbell PJ, Futreal PA. The cancer genome. Nature. Nature Publishing Group; 2009;458:719-24.

Storey JD, Tibshirani R. Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences of the United States of America. 2003;100:9440–5.
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    #Loading the Genomic Event data object derived from variants annotations
    #identified in 188 Lung Adenocarcinoma patients (Ding et al, 2008)
    data(LUAD_CaseStudy_updatedGS)
    
    #Loading KEGG pathway gene-set collection data object obtained by post-processing
    data(SLAPE.MSigDB_KEGG_hugoUpdated)
    
    #Loading genome-wide total exonic block lengths
    data(SLAPE.all_genes_exonic_content_block_lengths_ensemble)
    
    #Running SLAPenrich analysis
    PFPw<-SLAPE.analyse(EM = LUAD_CaseStudy_ugs,PATH_COLLECTION = KEGG_PATH,
                       show_progress = TRUE,
                       NSAMPLES = 1,
                       NGENES = 1,
                       accExLength = TRUE,
                       BACKGROUNDpopulation = rownames(LUAD_CaseStudy_ugs),
                       path_probability = 'Bernoulli',
                       GeneLenghts = GECOBLenghts)
    
    #Show the top-10 SLAPenriched pathway with an exclusive coverage > 80\%
    unlist(KEGG_PATH$PATHWAY[PFPw$pathway_id[which(PFPw$pathway_exclusiveCoverage>80)]][1:10])
  }
\keyword{Analysis}