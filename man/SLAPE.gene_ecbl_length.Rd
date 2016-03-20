\name{SLAPE.gene_ecbl_length}
\alias{SLAPE.gene_ecbl_length}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing the total exonic block length of a given gene
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SLAPE.gene_ecbl_length(ExonAttributes, GENE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExonAttributes}{Dataframe containing genomic coordinates of all the exon for all the genes in the genome. This is available in the
  \code{SLAPE.all_genes_exonic_lengths_ensemble_20160209} data object.
}
  \item{GENE}{
    A official HUGO gene symbol.
}
}
\details{

}
\value{
    An integer value specifying the total exonic block length of the genes specified in \code{GENE}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
All the genome-wide exonic block lengths are precomputed and available in the \code{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209} data object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.all_genes_exonic_lengths_ensemble_20160209}},
\code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
