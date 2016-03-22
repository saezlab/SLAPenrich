\name{SLAPE.heuristic_mut_ex_sorting}
\alias{SLAPE.heuristic_mut_ex_sorting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run-minimizing permutation of rows and colums of a binary matrix highlighting trends of mutual exclusivity (\emph{mutual exclusivity sorting})
}
\description{
The set of somatic mutations of a cancer genomic dataset can be easily modeled as a binary matrix, where columns indicate samples, and rows indicate genes (or vice-versa) and non-zero entries indicate the presence of somatic mutations in given gene/sample combinations. In a binary matrix, a run is a sequence of consecutive non-zero entries. Reordering the rows and the columns of a binary matrix modeling a sub-set of genes and samples of a genomic dataset in a way that the number of runs on its rows, as well as its column-wise marginal totals, are minimized is an effective way to highlight patterns of mutual exclusivity among the runs of different rows, therefore among the groups of samples harbouring mutations in at least one gene of the considered sub-set. Here we call such permutations of rows and columns a \emph{mutual exclusivity sorting} of the inputted binary matrix.
}
\usage{
SLAPE.heuristic_mut_ex_sorting(mutPatterns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutPatterns}{
    A generic binary matrix with row and column names.
}
}
\details{
This function implements a heruistic to perform the \emph{mutual exclusivity sorting} of a binary matrix in order to minimize the number row-wise runs and the column-wise marginal totals. To make the description of the algorithm simpler we will assume that the inputted binary matrix models a genomic dataset with samples on the columns and genes on the rows (as detailed in the description).

In the initial step of the algorithm all the samples and all the genes in the input matrix are declared \emph{uncovered} and an empty vector is initialized: this is the list of \emph{covered} genes \eqn{G}.

Then the algorithm proceeds through a series of iterations until the sets of uncovered genes and uncovered samples are both empty.

In each of these iterations a \emph{best in class gene} is identified. This is the uncovered gene with the maximal \emph{exclusive coverage}, which is defined as the number of uncovered samples in which this gene is mutated minus the number of samples in which at least another uncovered gene is mutated.

Finally, the identified best in class gene is removed from the set of the uncovered genes, it is attached to \eqn{G}, and the set of samples in which it is mutated are removed from the set of the uncovered samples.

After these iterations have been executed, an empty vector of samples \eqn{L} is initialized, all the samples of the dataset are labeled again as uncovered, an empty vector of samples \eqn{S} is initialized.

Then for each of the best in class gene \eqn{g} (in the same order as they appear in \eqn{G}) and until there are uncovered samples, the uncovered samples in which \eqn{g} is mutated are sorted according to the exclusive coverage of {g} across them (in decreasing ordered), they are labeled as covered samples and attached in the resulting order to \eqn{L}.

To obtain the final mutual-exclusivity sorting of the initial dataset, the corresponding inputted binary matrix is rearranged by permuting the genes/rows in the same order as they appear in \eqn{G} and the samples/columns in the same order as they appear in \eqn{L}.

}
\value{
    The inputted binary matrix with rows and columns permuted by the heuristic mutual exclusivity sorting algorithm described in the details.
}
\references{
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    #Generating a random sparse binary matrix
    dataset<-matrix(0,10,20,dimnames = list(paste('g',1:10,sep=''),paste('s',1:20,sep='')))
    dataset[sample(sample(200,30))]<-1

    #Visualising a heatmap of the original dataset,
    #with blue cells indicating non-null entries
    pheatmap(dataset,cluster_rows = FALSE,cluster_cols = FALSE,
         legend = FALSE,main='Original dataset',col=c('white','blue'))

    #wait 
    cat ("Press [enter] to continue")
    line <- readline()

    #Mutual exclusivity sorting the binary matrix
    me_sorted_dataset<-SLAPE.heuristic_mut_ex_sorting(dataset)

    #Visualising a heatmap of the mutual exclusivity sorted dataset,
    #with blue cells indicating non-null entries
    pheatmap(me_sorted_dataset,cluster_rows = FALSE,cluster_cols = FALSE,
             legend = FALSE,main='Original dataset',col=c('white','blue'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Result-exploration}

