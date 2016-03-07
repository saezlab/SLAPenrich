load('data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData')

bernoully_gene_mutProb<-function(gene,rho=10^-6){
    L<-GECOBLenghts[gene]
    b_i<-1-exp(-rho*L)
    return(b_i)
}

singleSampleExactProbMass<-function(gene_set,mutated_gene_set,rho=10^-6){
    
    b<-bernoully_gene_mutProb(mutated_gene_set,rho = rho)
    R<-(1-b)/b
    
    G<-sum(GECOBLenghts[mutated_gene_set])

    R<-prod(R)
    
    m<-length(gene_set)
    k<-length(mutated_gene_set)
    
    TR<-R
    
    for (Y in 1:k){
        upperLim<-(m-(Y-1))
        lowerLim<-k-Y+1
        TR<-TR*(upperLim-lowerLim+1)
    } 

    exp(-rho*G)*TR
}