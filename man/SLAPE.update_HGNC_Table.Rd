\name{SLAPE.update_HGNC_Table}
\alias{SLAPE.update_HGNC_Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Updating the R data object containing a HUGO approved catalogue of gene symbols and synomyms
}
\description{
This function creates a data frame containing up-to-date HUGO Gene Nomenclature Committee (HGNC) approved symbols and their synonyms downloading updated relevant information from the HUGO Gene Nomenclature Committee web-portal (\url{http://www.genenames.org}).
This table is used by the
\code{\link{SLAPE.check_and_fix_gs_Dataset}} and \code{\link{SLAPE.check_and_fix_path_collection}} functions to check and update the gene symbol identifiers of a integer genomic event matrix. A precomputed data frame (created on February 2016) is available in the
\code{\link{SLAPE.hgnc.table}} data object.
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
\code{\link{SLAPE.hgnc.table}}

\code{\link{SLAPE.check_and_fix_path_collection}}

\code{\link{SLAPE.check_and_fix_gs_Dataset}}

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
