\name{SLAPE.update_HGNC_Table}
\alias{SLAPE.update_HGNC_Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SLAPE.update_HGNC_Table
}
\description{
This function creates a data frame containing up-to-date HUGO Gene Nomenclature Committee (HGNC) approved symbols and their synonyms downloading updated relevant information from the HUGO Gene Nomenclature Committee web-portal (http://www.genenames.org).
This table is used by the \code{SLAPE.update_HGNC_Table} function to check and update the gene symbol identifiers of a integer genomic event matrix. A precomputed data frame is available in the \code{SLAPE.hgnc.table_20160210} data object.
}
\usage{
SLAPE.update_HGNC_Table()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    print("Downloading updated table from genaname.org")
    day <- Sys.time()
    day <- str_split(day, " ")[[1]][1]
    fn <- paste("externalData/hgnc_complete_set_", day, ".txt", 
        sep = "")
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt", 
        destfile = fn)
    fc <- file(fn)
    mylist <- str_split(readLines(fc), "\t")
    close(fc)
    mylist <- mylist[2:length(mylist)]
    ngenes <- length(mylist)
    flag <- TRUE
    print("Updating HGNC table...")
    pb <- txtProgressBar(min = 0, max = ngenes, initial = 0, 
        style = 3)
    for (i in 1:ngenes) {
        setTxtProgressBar(pb, i)
        if (length(grep("withdrawn", mylist[[i]])) == 0) {
            currentApprovedSymbol <- mylist[[i]][2]
            currentAlias <- mylist[[i]][9]
            currentAlias <- setdiff(unlist(str_split(currentAlias, 
                "\"")), "")
            currentAlias <- unlist(str_split(currentAlias, "[|]"))
            currentPrevSymbol <- mylist[[i]][11]
            currentPrevSymbol <- setdiff(unlist(str_split(currentPrevSymbol, 
                "\"")), "")
            currentPrevSymbol <- unlist(str_split(currentPrevSymbol, 
                "[|]"))
            prevSymbAndAliases <- c(currentApprovedSymbol, currentAlias, 
                currentPrevSymbol)
            currentChunk <- cbind(sort(prevSymbAndAliases), rep(currentApprovedSymbol, 
                length(prevSymbAndAliases)))
            if (flag) {
                updated.hgnc.table <- currentChunk
                flag <- FALSE
            }
            else {
                updated.hgnc.table <- rbind(updated.hgnc.table, 
                  currentChunk)
            }
        }
    }
    close(pb)
    print("DONE!")
    colnames(updated.hgnc.table) <- c("Symbol", "Approved.Symbol")
    rownames(updated.hgnc.table) <- as.character(1:nrow(updated.hgnc.table))
    updated.hgnc.table <- data.frame(updated.hgnc.table)
    return(updated.hgnc.table)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data-management}
