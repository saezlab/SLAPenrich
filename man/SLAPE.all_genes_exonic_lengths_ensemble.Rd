\name{SLAPE.all_genes_exonic_lengths_ensemble}
\alias{SLAPE.all_genes_exonic_lengths_ensemble}
\docType{data}
\title{
Genome-wide exone attributes and genomic coordinates
}
\description{
A data frame containing attributes and chromosomal coordinate of all the gene exones
}
\format{
  A data frame with 553609 rows (one for each exone) and the following columns
  \describe{
    \item{\code{ensembl_gene_id}}{String vector containing ensemble gene identifiers;}
    \item{\code{external_gene_name}}{String vector containing gene names;}
    \item{\code{exon_chrom_start}}{Numerical vector containing chromosomal start positions;}
    \item{\code{exon_chrom_end}}{Numerical vector containing chromosomal end positions}
  }
  Variable name: \code{GEA}.
}
\details{
This list has been assembled by using functions from the \code{biomaRt} package and can be updated 
using the \code{SLAPE.update_exon_attributes}.
}
\note{This object is used by the \code{SLAPE.compute_gene_exon_content_block_lengths} and \code{SLAPE.gene_ecbl_length} function to compute
the genome-wide total exonic block lengths or the total exonic block length of a given gene, respectively.}
\seealso{
\code{\link{SLAPE.update_exon_attributes}},
\code{\link{SLAPE.compute_gene_exon_content_block_lengths}}
\code{\link{SLAPE.gene_ecbl_length}}
}
\keyword{datasets}
