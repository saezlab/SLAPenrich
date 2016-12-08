\name{SLAPE.hgnc.table}
\alias{SLAPE.hgnc.table}
\docType{data}
\title{
HUGO gene symbols and their previous synonyms (up to February 2016)
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("SLAPE.hgnc.table")}
\format{
  A data frame with approved HUGO gene symbols in one column \code{Approved.Symbol} and their previously approved synonyms \code{Symbol} in another column (up to February 2016).
  Variable Name: \code{hgnc.table}.
}
\source{
HUGO Gene Nomenclature Committee web-portal (\url{http://www.genenames.org})
}
\note{
    This object can be updated to a more recent version using the \code{SLAPE.update_HGNC_Table} function.
}

\seealso{
\code{\link{SLAPE.update_HGNC_Table}}
}
\examples{
data(SLAPE.hgnc.table)
## maybe str(SLAPE.hgnc.table_20160210) ; plot(SLAPE.hgnc.table_20160210) ...
}
\keyword{datasets}
