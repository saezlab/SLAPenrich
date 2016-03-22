\name{SLAPE.update_HGNC_Table}
\alias{SLAPE.update_HGNC_Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Updating the catalogue of approved gene symbols and synomyms
}
\description{
This function creates a data frame containing up-to-date HUGO Gene Nomenclature Committee (HGNC) approved symbols and their synonyms downloading updated relevant information from the HUGO Gene Nomenclature Committee web-portal (\url{http://www.genenames.org}).
This table is used by the
\code{\link{SLAPE.check_and_fix_gs_Dataset}} function to check and update the gene symbol identifiers of a integer genomic event matrix. A precomputed data frame is available in the
\code{\link{SLAPE.hgnc.table_20160210}} data object.
}
\usage{
SLAPE.update_HGNC_Table()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with updated approved HUGO gene symbols in one column and their previously approved synonyms in another column.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
