\name{SLAPE.hgnc.table_20160210}
\alias{SLAPE.hgnc.table_20160210}
\docType{data}
\title{
Recent HUGO gene symbols and their previous synonyms
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("SLAPE.hgnc.table_20160210")}
\format{
  A data frame with updated approved HUGO gene symbols in one column \code{Approved.Symbol} and their previously approved synonyms \code{Symbol} in another column.
  Please note that the name of this data frame is \code{updated.hgnc.table}.
}
\source{
HUGO Gene Nomenclature Committee web-portal (\url{http://www.genenames.org})
}
\note{
    This object can be updated using the \code{SLAPE.update_HGNC_Table} function.
}

\seealso{
\code{\link{SLAPE.update_HGNC_Table}}
}
\examples{
data(SLAPE.hgnc.table_20160210)
## maybe str(SLAPE.hgnc.table_20160210) ; plot(SLAPE.hgnc.table_20160210) ...
}
\keyword{datasets}
