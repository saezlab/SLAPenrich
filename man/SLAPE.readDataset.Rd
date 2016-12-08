\name{SLAPE.readDataset}
\alias{SLAPE.readDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading genomic evant dataset
}
\description{
This function reads a genomic dataset from a csv file and it stores it into an integer matrix.
Row names of this matrix are official gene symbols and column names are sample identifiers.
A non-zero entry in the \eqn{i,j} position indicates the presence of  somatic mutation hosted by the \eqn{i}-th gene in the \eqn{j}-th sample (if the matrix is binary) or the number of point mutations hosted by the \eqn{i}-th gene in the \eqn{j}-th sample (if the matrix contains integers).
}
\usage{
SLAPE.readDataset(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The path of the csv file to be read and stored in the genomic event matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An integer matrix modeling a genomic event dataset. Row names of this matrix are official gene symbols and column names are sample identifiers.
A non-zero entry in the \eqn{i,j} position indicates the presence of a somatic mutation hosted by the \eqn{i}-th gene in the \eqn{j}-th sample
(if the matrix is binary) or the number of point mutations hosted by the \eqn{i}-th gene in the \eqn{j}-th sample (if the matrix contains integers).
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
