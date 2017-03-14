\name{SLAPE.check_and_fix_path_collection}
\alias{SLAPE.check_and_fix_path_collection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check and fix gene symbol names in a collection of pathway gene sets.
}
\description{
This function checks that gene identifiers contained in a pathway gene set collection are actually updated official gene symbols approved by the HUGO Gene Nomenclature Committee (HGNC)
(http://www.genenames.org).
}
\usage{
SLAPE.check_and_fix_path_collection(pathColl, updated.hgnc.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathColl}{
  A list containing pathway gene-sets and annotations.
}
  \item{updated.hgnc.table}{
    A data frame containing up-to-date approved HGNC symbols (Approved.Symbol variable) and their synonyms (Symbol variable). This is availalbe in the
    \code{\link{SLAPE.hgnc.table}}
    data object or can be created by downloading updated relevant information from the HUGO Gene Nomenclature Committee web-portal (http://www.genenames.org), using the function \code{\link{SLAPE.update_HGNC_Table}}.
    
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Pathway collection provided in input but with gene identifiers updated to the most recent approved gene symbols and
    not approved symbols removed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Iorio - iorio@ebi.ac.uk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SLAPE.update_HGNC_Table}},
\code{\link{SLAPE.hgnc.table}}
}
\examples{
data(SLAPE.PATHCOM_HUMAN)
data(SLAPE.hgnc.table)
updatedGeneSymbolsDataset<-
SLAPE.check_and_fix_path_collection(PATHCOM_HUMAN,hgnc.table)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
