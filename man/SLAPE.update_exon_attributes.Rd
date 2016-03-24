\name{SLAPE.update_exon_attributes}
\alias{SLAPE.update_exon_attributes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating an updated gene exon attributes data object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SLAPE.update_exon_attributes()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe containing updated genomic coordinates of all the exon for all the genes in the genome, from Ensemble.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
A dataframe containing genomic coordinates of all the exon for all the genes in the genome, from Ensemble is precomputed and available in the
  \code{\link{SLAPE.all_genes_exonic_lengths_ensemble_20160209}} data object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209}},
\code{\link{SLAPE.all_genes_exonic_lengths_ensemble_20160209}},
\code{\link{SLAPE.compute_gene_exon_content_block_lengths}},
\code{\link{SLAPE.gene_ecbl_length}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}

