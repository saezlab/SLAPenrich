\name{SLAPE.compute_gene_exon_content_block_lengths}
\alias{SLAPE.compute_gene_exon_content_block_lengths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing genome-wide total exonic block lengths
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SLAPE.compute_gene_exon_content_block_lengths(ExonAttributes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExonAttributes}{
    Dataframe containing genomic coordinates of all the exon for all the genes in the genome. This is available in the
    \code{\link{SLAPE.all_genes_exonic_lengths_ensemble_20160209}} data object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A genome-wide named vector containing the total exonic block lengths for all the genes. Names of this vector are offical HUGO gene symbols.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
The genome-wide exonic block lengths are precomputed and available in the
\code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209}}
data object, this function can be used to update this data object with the most-up-to-date information from Ensemble (using biomaRt functions).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209}},
\code{\link{SLAPE.all_genes_exonic_lengths_ensemble_20160209}},
\code{\link{SLAPE.compute_gene_exon_content_block_lengths}},
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}

