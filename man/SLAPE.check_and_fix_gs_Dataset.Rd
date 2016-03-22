\name{SLAPE.check_and_fix_gs_Dataset}
\alias{SLAPE.check_and_fix_gs_Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check and fixing genomic event dataset
}
\description{
This function checks that the row names of a genomic dataset are actually updated official gene symbols approved by the HUGO Gene Nomenclature Committee (HGNC)
(http://www.genenames.org).
}
\usage{
SLAPE.check_and_fix_gs_Dataset(Dataset, updated.hgnc.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dataset}{
  An integer matrix modeling a genomic event dataset where row names are gene symbols and column names are sample identifiers.
A non-null entry in the i,j position indicates the presence of a somatic mutation hosted by the i-th gene in the j-th sample
(if the matrix is binary) or the number of point mutations hosted by the i-th gene in the j-th sample (if the matrix contains integers).
}
  \item{updated.hgnc.table}{
    A data frame containing up-to-date approved HGNC symbols (Approved.Symbol variable) and their synonyms (Symbol variable). This is availalbe in the
    \code{\link{SLAPE.hgnc.table_20160210}}
    data object or can be created by downloading updated relevant information from the HUGO Gene Nomenclature Committee web-portal (http://www.genenames.org), using the function \code{\link{SLAPE.update_HGNC_Table}}.
    
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    The integer matrix provided in input but with row names updated to the most recent approved gene symbol and rows with not found
    gene synonyms as names removed.
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
\code{\link{SLAPE.hgnc.table_20160210}}
}
\examples{
data(LUAD_CaseStudy)
data(SLAPE.hgnc.table_20160210)
updatedGeneSymbolsDataset<-SLAPE.check_and_fix_gs_Dataset(LUAD_CaseStudy,updated.hgnc.table)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
