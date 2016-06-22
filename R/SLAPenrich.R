# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#internal functions
SLE.hypTest<-function(x,k,n,N){
    PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
    return(PVALS)
}
SLE.rearrangeMatrix<-function(patterns,GENES){
    
    remainingSamples<-colnames(patterns)
    
    toAdd<-NULL
    
    for (g in GENES){
        remainingGenes<-setdiff(GENES,g)
        
        P1<-matrix(c(patterns[g,remainingSamples]),length(g),length(remainingSamples),dimnames = list(g,remainingSamples))
        P2<-matrix(c(patterns[remainingGenes,remainingSamples]),length(remainingGenes),length(remainingSamples),
                   dimnames=list(remainingGenes,remainingSamples))
        
        if(length(remainingGenes)>1){
            DD<-colnames(P1)[order(P1-colSums(P2),decreasing=TRUE)]
        }else{
            DD<-colnames(P1)[order(P1-P2,decreasing=TRUE)]
        }
        
        toAdd<-c(toAdd,names(which(patterns[g,DD]>0)))
        remainingSamples<-setdiff(remainingSamples,toAdd)
        if(length(remainingSamples)==0){
            break
        }
    }
    
    toAdd<-c(toAdd,remainingSamples)
    
    return(toAdd)
}
SLE.findBestInClass<-function(patterns){
    
    if(nrow(patterns)==1){
        return(rownames(patterns))
    }
    
    if(ncol(patterns)==1){
        return(rownames(patterns)[1])
    }
    
    genes<-rownames(patterns)
    
    exclCov<-rep(NA,length(genes))
    names(exclCov)<-genes
    for (g in genes){
        residGenes<-setdiff(genes,g)
        if (length(residGenes)>1){
            exclCov[g]<-sum(patterns[g,]-colSums(patterns[residGenes,]))
        }else{
            exclCov[g]<-sum(patterns[g,]-patterns[residGenes,])
        }
    }
    
    return(names(sort(exclCov,decreasing=TRUE))[1])
}
SLE.buildPathMemb<-function(GENES,pathway_ids,PATH_COLLECTION){
    np<-length(pathway_ids)
    ngenes<-length(GENES)
    MM<-matrix(0,ngenes,np,dimnames = list(GENES,pathway_ids))
    
    for (i in 1:np){
        MM[,i]<-is.element(GENES,PATH_COLLECTION$HGNC_SYMBOL[[pathway_ids[i]]])+0
    }
    
    return(MM)
}
SLE.assignTotalLengthTOpathways<-function(){
    np<-length(PATH_COLLECTION$PATHWAY)
    
    Glenghts<-list()
    TOTAL_G_LENGHT<-vector()
    for (i in 1:np){
        print(i)
        Glenghts[[i]]<-gene_in_path_lengths[PATH_COLLECTION$HGNC_SYMBOL[[i]]]
        TOTAL_G_LENGHT[i]<-sum(gene_in_path_lengths[PATH_COLLECTION$HGNC_SYMBOL[[i]]],na.rm = TRUE)
    }
    
    names(Glenghts)<-PATH_COLLECTION$PATHWAY
    names(TOTAL_G_LENGHT)<-PATH_COLLECTION$PATHWAY
    
    PATH_COLLECTION$Glengths<-Glenghts
    PATH_COLLECTION$TOTAL_G_LENGHT<-TOTAL_G_LENGHT
    
    return(PATH_COLLECTION)
    
}
SLE.combine<-function(wResults,unwResults,fdrTH=20){
    
    u_ids<-unwResults$pathway_id[which(unwResults$pathway_perc_fdr<=fdrTH)]
    w_ids<-wResults$pathway_id[which(wResults$pathway_perc_fdr<=fdrTH)]
    
    common_ids<-intersect(u_ids,w_ids)
    
    cResults<-list()
    
    cResults$pathway_id<-common_ids
    
    
    u_ids<-match(common_ids,unwResults$pathway_id)
    w_ids<-match(common_ids,wResults$pathway_id)
    
    
    cResults$W_pathway_perc_fdr<-wResults$pathway_perc_fdr[w_ids]
    cResults$U_pathway_perc_fdr<-unwResults$pathway_perc_fdr[u_ids]
    
    tmp<-order(rowMeans(cbind(cResults$U_pathway_perc_fdr,cResults$W_pathway_perc_fdr)))
    
    cResults$W_pathway_perc_fdr<-wResults$pathway_perc_fdr[w_ids[tmp]]
    cResults$U_pathway_perc_fdr<-unwResults$pathway_perc_fdr[u_ids[tmp]]
    
    cResults$U_pathway_logOddRatios<-unwResults$pathway_logOddRatios[u_ids[tmp]]
    cResults$W_pathway_logOddRatios<-wResults$pathway_logOddRatios[w_ids[tmp]]
    
    cResults$U_pathway_pvals<-unwResults$pathway_pvals[u_ids[tmp]]
    cResults$W_pathway_pvals<-wResults$pathway_pvals[w_ids[tmp]]
    
    return(cResults)
}
SLE.plotMyHeat <- function(x,orPlot,verdata,filename) {
    pdf(filename,height = 15,width=15)
    l <- matrix(c(1,2), ncol = 2)
    layout(l)
    op <- par(mar = c(30,2,1,0))
    
    if(ncol(x)>1){COLS<-c('white','darkgreen')}
    else{COLS<-c('darkgreen')}
    image(t(x),axes=FALSE,col=COLS)
    
    NAMES<-names(verdata)
    NAMES<-str_sub(NAMES,1,50)
    
    stverdata<-rep('',length(verdata))
    for (j in 1:length(verdata)){
        if (verdata[j]<1){
            stverdata[j]<-format(verdata[j],scientific=TRUE,digits=2)
        }else{
            stverdata[j]<-format(verdata[j],scientific=FALSE,digits=2)
        }
    }
    
    NAMES<-paste(NAMES,' (FDR ',stverdata,'%)',sep='')
    
    if(ncol(x)>1){
        axis(side = 1,at = seq(0,1,1/(ncol(x)-1)),las=2,labels=NAMES)
    }else{
        axis(side = 1,at = 1,las=2,labels=NAMES)
    }
    
    
    par(mar = c(30,6,1,1))
    barplot(orPlot,horiz = TRUE, yaxs = "i",las=2,xlim=c(0,max(orPlot)+10),cex.names = 0.8)
    
    par(op)
    dev.off()
}


#exported functions
SLAPE.check_and_fix_path_collection<-function(pathColl,updated.hgnc.table){
    
    
    npath<-length(pathColl$PATHWAY)
    
    for (i in 1:npath){
        print(i)
        tmp<-pathColl$HGNC_SYMBOL[[i]]
        
        if (length(tmp)>0){
            
        
        checked_gs<-checkGeneSymbols(tmp,hgnc.table = updated.hgnc.table)
        non_approved_id<-which(checked_gs[,2]==FALSE)
        
        if (length(non_approved_id)>0){
            
            print('The dataset contains non-approved gene symbols...')
            print('Outdated gene symbols have been updated:')
            
            outdated_id<-non_approved_id[which(!is.na(checked_gs[non_approved_id,3]))]
            print(paste(checked_gs[outdated_id,1],'->',checked_gs[outdated_id,3]))
            
            non_approved_id<-which(is.na(checked_gs[,3]))
            if(length(non_approved_id)>0){
                print('The following non approved gene symbols have been removed:')
                print(checked_gs[non_approved_id,1])
                
            }
            
            tmp[outdated_id]<-checked_gs[outdated_id,3]
            tmp<-setdiff(tmp,tmp[non_approved_id])
            
            pathColl$HGNC_SYMBOL[[i]]<-tmp
        }
            pathColl$Ngenes[[i]]<-length(pathColl$HGNC_SYMBOL[[i]])
        }
    }
    
    ui<-unlist(pathColl$Ngenes)
    idx<-which(ui>2)
    
    pathColl$PATHWAY<-pathColl$PATHWAY[idx]
    pathColl$HGNC_SYMBOL<-pathColl$HGNC_SYMBOL[idx]
    pathColl$Hallmark<-pathColl$Hallmark[idx]
    pathColl$Ngenes<-pathColl$Ngenes[idx]
    pathColl$backGround<-sort(unique(unlist(pathColl$HGNC_SYMBOL)))
    
    return(pathColl)
}
#documented
SLAPE.readDataset<-function(filename){
    fc<-as.matrix(read.csv(filename,row.names=1))
    fc[is.na(fc)]<-0
    
    return(fc)
}
SLAPE.check_and_fix_gs_Dataset<-function(Dataset,updated.hgnc.table){
    checked_gs<-checkGeneSymbols(rownames(Dataset),hgnc.table = updated.hgnc.table)    

    non_approved_id<-which(checked_gs[,2]==FALSE)
    
    if (length(non_approved_id)>0){

        print('The dataset contains non-approved gene symbols...')
        print('Outdated gene symbols have been updated:')

        outdated_id<-non_approved_id[which(!is.na(checked_gs[non_approved_id,3]))]
        print(paste(checked_gs[outdated_id,1],'->',checked_gs[outdated_id,3]))
        
        non_approved_id<-which(is.na(checked_gs[,3]))
        if(length(non_approved_id)>0){
            print('The following non approved gene symbols have been removed:')
            print(checked_gs[non_approved_id,1])
        
            }
        rownames(Dataset)[outdated_id]<-checked_gs[outdated_id,3]
        Dataset<-Dataset[setdiff(rownames(Dataset),rownames(Dataset)[non_approved_id]),]
        }
    
    return(Dataset)
}
SLAPE.update_HGNC_Table<-function(){
    print('Downloading updated table from genaname.org')
    day<-Sys.time()
    day<-str_split(day,' ')[[1]][1]
    fn<-paste('externalData/hgnc_complete_set_',day,'.txt',sep='')
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt",
                  destfile = fn)
    
    fc <- file(fn)
    mylist <- str_split(readLines(fc), "\t")
    close(fc)
    
    mylist<-mylist[2:length(mylist)]
    
    ngenes<-length(mylist)
    
    flag<-TRUE
    print('Updating HGNC table...')
    pb<-txtProgressBar(min = 0, max = ngenes, initial = 0,style = 3)
    
    for (i in 1:ngenes){
        setTxtProgressBar(pb,i)
        if(length(grep('withdrawn',mylist[[i]]))==0){
            currentApprovedSymbol<-mylist[[i]][2]
            currentAlias<-mylist[[i]][9]
            currentAlias<-setdiff(unlist(str_split(currentAlias,'"')),'')
            currentAlias<-unlist(str_split(currentAlias,'[|]'))
            
            currentPrevSymbol<-mylist[[i]][11]
            currentPrevSymbol<-setdiff(unlist(str_split(currentPrevSymbol,'"')),'')
            currentPrevSymbol<-unlist(str_split(currentPrevSymbol,'[|]'))
            
            prevSymbAndAliases<-c(currentApprovedSymbol,currentAlias,currentPrevSymbol)
            
            currentChunk<-cbind(sort(prevSymbAndAliases),rep(currentApprovedSymbol,length(prevSymbAndAliases)))
            
            if(flag){
                updated.hgnc.table<-currentChunk
                flag<-FALSE
            }else{
                updated.hgnc.table<-rbind(updated.hgnc.table,currentChunk)
            }
        }
        
    }
    close(pb)
    print('DONE!')
    
    colnames(updated.hgnc.table)<-c("Symbol","Approved.Symbol")
    rownames(updated.hgnc.table)<-as.character(1:nrow(updated.hgnc.table))
    
    
    updated.hgnc.table<-data.frame(updated.hgnc.table)
    return(updated.hgnc.table)
}
SLAPE.gene_ecbl_length<-function(ExonAttributes,GENE){
    
    id<-which(is.element(ExonAttributes$external_gene_name,GENE))
    
    if (length(id)==1){
        ECBlenght<-ExonAttributes$exon_chrom_end[id]-ExonAttributes$exon_chrom_start[id]+1
    }else{
        startPos<-ExonAttributes$exon_chrom_start[id]
        endPos<-ExonAttributes$exon_chrom_end[id]
        
        minPos<-min(startPos)
        maxPos<-max(endPos)
        
        span<-maxPos-minPos+1
        
        nsegments<-length(startPos)
        
        cocMat<-matrix(0,nsegments,span,dimnames = list(1:nsegments,as.character(minPos:maxPos)))
        
        for (i in 1:nsegments){
            cocMat[i,as.character(startPos[i]:endPos[i])]<-1
        }
        
        ECBlenght<-sum(sign(colSums(cocMat)))
        
    }
    
    return(ECBlenght)
}
SLAPE.compute_gene_exon_content_block_lengths<-function(ExonAttributes){
    print("Updating Genome-Wide Exon content block lengths... please wait...")
    uniqueGS<-unique(ExonAttributes$external_gene_name)
    ngenes<-length(uniqueGS)
    
    GECOBLenghts<-rep(NA,ngenes)
    names(GECOBLenghts)<-uniqueGS[1:ngenes]
    
    pb<-txtProgressBar(min = 0, max = ngenes, initial = 0,style = 3)
    
    for (i in 1:ngenes){
        setTxtProgressBar(pb,i)
        GECOBLenghts[i]<-SLAPE.gene_ecbl_length(ExonAttributes,uniqueGS[i])
    }
    close(pb)
    print("DONE!")
    return(GECOBLenghts)
}
SLAPE.update_exon_attributes<-function(){
    print("Updating Gene Exons Attributes... please wait...")
    EnsemblMart = useEnsembl(biomart="ensembl",dataset = 'hsapiens_gene_ensembl')
    
    
    print('     - Genome-wide list of HGNC coded genes: Download in progress...')
    HGNC_coded_Genes<-sort(unique(getBM(mart=EnsemblMart,attributes = c("ensembl_gene_id",
                                                                        "external_gene_name"),
                                        filters="with_hgnc",values = TRUE)[,2]))
    print('     + DONE!')
    
    ngenes<-length(HGNC_coded_Genes)
    nblocks<-floor(ngenes/100)
    
    print('     - Downloading exon attributes...')
    
    pb<-txtProgressBar(min = 0, max = nblocks+1, initial = 0,style = 3)
    
    for (i in 1:nblocks){
        setTxtProgressBar(pb,i)
        currentExons<-getBM(mart=EnsemblMart,attributes = c("ensembl_gene_id",
                                                            "external_gene_name",
                                                            "exon_chrom_start",
                                                            "exon_chrom_end"),
                            filters="hgnc_symbol",values = HGNC_coded_Genes[(100*(i-1)+1):(100*i)])
        
        if (i==1){
            ExonAttributes<-currentExons
        }else{
            ExonAttributes<-rbind(ExonAttributes,currentExons)   
        }
    }
    
    remainingGenes<-ngenes-nblocks*100
    currentExons<-getBM(mart=EnsemblMart,attributes = c("ensembl_gene_id",
                                                        "external_gene_name",
                                                        "exon_chrom_start",
                                                        "exon_chrom_end"),
                        filters="hgnc_symbol",values = HGNC_coded_Genes[(100*i+1):(100*i+remainingGenes)])
    
    ExonAttributes<-rbind(ExonAttributes,currentExons)   
    setTxtProgressBar(pb,i+1)
    close(pb)
    print('     - DONE!')
    
    print("DONE!")
    return(ExonAttributes)
}
SLAPE.analyse<-function(EM,
                        show_progress=TRUE,
                        correctionMethod='fdr',
                        NSAMPLES=1,
                        NGENES=1,
                        accExLength=TRUE,
                        BACKGROUNDpopulation=NULL,
                        PATH_COLLECTION,
                        path_probability='Bernoulli',
                        GeneLenghts,RHO=NA){
    
    if(length(BACKGROUNDpopulation)>0){
        if(length(which(duplicated(BACKGROUNDpopulation)))>0){
            warning('Duplicated gene names found in the background population')
        }
        BACKGROUNDpopulation<-unique(BACKGROUNDpopulation)
    }
    
    if(sum(is.na(EM))>0){
        error('The inputted EM contains NA!')
    }
    
    if(accExLength){
        cat('Sample fingerprinting for pathawy alterations (accounting for gene exonic length)...\n')
        gnames<-unlist(PATH_COLLECTION$HGNC_SYMBOL)
        if(!length(BACKGROUNDpopulation)){
            N<-sum(GeneLenghts[unique(gnames)],na.rm = TRUE)
        }else{
            N<-sum(GeneLenghts[unique(BACKGROUNDpopulation)],na.rm = TRUE)
        }
    }else{
        cat('Sample fingerprinting for pathawy alterations...\n')
        if(!length(BACKGROUNDpopulation)){
            N<-unlist(PATH_COLLECTION$HGNC_SYMBOL)
            N<-length(unique(N))
        }else{
            N<-length(BACKGROUNDpopulation)
        }
    }
    
    np<-length(PATH_COLLECTION$PATHWAY)
    
    nsamples<-ncol(EM)
    
    PN<-PATH_COLLECTION$PATHWAY
    
    pathway_EM<-matrix(0,nrow=np,ncol=nsamples,dimnames=list(1:np,colnames(EM)))
    pathway_Probability<-matrix(0,nrow=np,ncol=nsamples,dimnames=list(1:np,colnames(EM)))
    pathway_grandTotal<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_pvals<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_mus<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathwayExclusive_coverage<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_cometPvalue<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_cometFDR<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_individualEMs<-list()
    
    if(show_progress){
        pb <- txtProgressBar(min=1,max=nsamples*np,style=3)
    }
    
    toExclude<-rep(FALSE,np)
    for (i in 1:np){
        
        currentGeneSet<-intersect(PATH_COLLECTION$HGNC_SYMBOL[[i]],rownames(EM))
        
        if(accExLength){
            GLENGHTS<-
                GeneLenghts[currentGeneSet]
            
            if (length(which(!is.na(GLENGHTS)))==0 | length(currentGeneSet)<=NGENES){
                toExclude[i]<-TRUE
            }
            n <- sum(GeneLenghts[currentGeneSet],na.rm = TRUE)
        }else{
            n <-length(currentGeneSet)
        }
        
        
        if (length(currentGeneSet)>1){
            pathway_individualEMs[[i]]<-EM[currentGeneSet,]
        }else{
            pathway_individualEMs[[i]]<-matrix(EM[currentGeneSet,],1,ncol(EM),dimnames = list(currentGeneSet,colnames(EM)))
        }
        
        sampleCovered<-sum(colSums(pathway_individualEMs[[i]])>0)
        exclusiveSampleCovered<-sum(colSums(pathway_individualEMs[[i]])==1)
        
        pathwayExclusive_coverage[i]<-100*exclusiveSampleCovered/sampleCovered
        
        if (length(currentGeneSet)>0){
            pathway_grandTotal[i]<-sum(unlist(c(EM[currentGeneSet,])))
        }else{
            pathway_grandTotal[i]<-0
        }
        
        for (j in 1:nsamples){
            k <- sum(EM[,j])
            x<-1
            if(length(currentGeneSet)>0){
                pathway_EM[i,j]<-sign(sum(EM[currentGeneSet,j]))
            }
            
            rho<-RHO
            
            if (is.na(rho)){
                rho<-k/N    
            }
            
            
            if (path_probability=='Bernoulli'){
                pathway_Probability[i,j]<-1-exp(-rho*n)    
            }
            if (path_probability=='HyperGeom'){
                pathway_Probability[i,j]<-SLE.hypTest(x,k,n,N)
            }
            
            if(show_progress){setTxtProgressBar(pb, (i-1)*nsamples+j)}
        }
    }
    
    if(show_progress){
        Sys.sleep(1)
        close(pb)
    }
    
    pathway_id<-1:length(PATH_COLLECTION$PATHWAY)
    
    pathway_mus<-rowSums(pathway_Probability)
    pathway_var<-rowSums(1-pathway_Probability)*pathway_Probability
    
    pathway_logOddRatios<-log10(rowSums(pathway_EM)/pathway_mus)
    
    for (i in 1:nrow(pathway_Probability)){
        pathway_pvals[i]<-sum(dpoibin(rowSums(pathway_EM)[i]:ncol(EM),pathway_Probability[i,]))
    }
    
    id<-which(rowSums(pathway_EM)>=NSAMPLES)
    
    pathway_id<-pathway_id[id]
    pathway_EM<-pathway_EM[id,]
    pathway_Probability<-pathway_Probability[id,]
    pathway_pvals<-pathway_pvals[id]
    pathway_mus<-pathway_mus[id]
    pathway_logOddRatios<-pathway_logOddRatios[id]
    pathway_individualEMs<-pathway_individualEMs[id]
    pathwayExclusive_coverage<-pathwayExclusive_coverage[id]
    toExclude<-toExclude[id]
    
    pathway_id<-pathway_id[order(pathway_pvals)]
    pathway_EM<-pathway_EM[order(pathway_pvals),]
    pathway_Probability<-pathway_Probability[order(pathway_pvals),]
    pathway_logOddRatios<-pathway_logOddRatios[order(pathway_pvals)]
    pathway_mus<-pathway_mus[order(pathway_pvals)]
    pathway_individualEMs<-pathway_individualEMs[order(pathway_pvals)]
    pathwayExclusive_coverage<-pathwayExclusive_coverage[order(pathway_pvals)]
    toExclude<-toExclude[order(pathway_pvals)]
    pathway_pvals<-sort(pathway_pvals)
    
    pathway_id<-pathway_id[which(!toExclude)]
    pathway_EM<-pathway_EM[which(!toExclude),]
    pathway_Probability<-pathway_Probability[which(!toExclude),]
    pathway_logOddRatios<-pathway_logOddRatios[which(!toExclude)]
    pathway_mus<-pathway_mus[which(!toExclude)]
    pathway_individualEMs<-pathway_individualEMs[which(!toExclude)]
    pathwayExclusive_coverage<-pathwayExclusive_coverage[which(!toExclude)]
    pathway_pvals<-pathway_pvals[which(!toExclude)]
    toExclude<-toExclude[which(!toExclude)]
    
    if (correctionMethod=='qvalue'){
        Q<-qvalue(pathway_pvals)
        pathway_perc_fdr<-100*Q$qvalues
    }else{
        pathway_perc_fdr<-100*p.adjust(pathway_pvals,method = correctionMethod)
    }
    
    cat('\nDone!...\n')
    
    names(pathway_pvals)<-names(pathway_logOddRatios)
    names(pathway_perc_fdr)<-names(pathway_logOddRatios)
    names(pathwayExclusive_coverage)<-names(pathway_logOddRatios)
    names(pathway_individualEMs)<-names(pathway_logOddRatios)
    
    return(list(pathway_id=pathway_id,
                pathway_EM=pathway_EM,
                pathway_Probability=pathway_Probability,
                pathway_mus=pathway_mus,
                pathway_logOddRatios=pathway_logOddRatios,
                pathway_pvals=pathway_pvals,
                pathway_perc_fdr=pathway_perc_fdr,
                pathway_exclusiveCoverage=pathwayExclusive_coverage,
                pathway_individualEMs=pathway_individualEMs))
}

SLAPE.write.table<-function(PFP,EM,filename='',fdrth=Inf,exclcovth=0,PATH_COLLECTION,GeneLenghts){
    
    id<-which(PFP$pathway_exclusiveCoverage>exclcovth & PFP$pathway_perc_fdr<fdrth)
    
    PFP$pathway_id<-PFP$pathway_id[id]
    PFP$pathway_mus<-PFP$pathway_mus[id]
    PFP$pathway_EM<-PFP$pathway_EM[id,]
    PFP$pathway_logOddRatios<-PFP$pathway_logOddRatios[id]
    PFP$pathway_pvals<-PFP$pathway_pvals[id]
    PFP$pathway_perc_fdr<-PFP$pathway_perc_fdr[id]
    PFP$pathway_exclusiveCoverage<-PFP$pathway_exclusiveCoverage[id]
    PFP$pathway_Probability<-PFP$pathway_Probability[id,]
    
    NAMES<-PATH_COLLECTION$PATHWAY[PFP$pathway_id]
    
    NAMES<-str_replace_all(NAMES,',','//')
    
    
    TOTexlength<-rep(NA,length(NAMES))
    
    for (i in 1:length(NAMES)){
        TOTexlength[i]<-sum(GeneLenghts[PATH_COLLECTION$HGNC_SYMBOL[[PFP$pathway_id[i]]]],na.rm=TRUE)    
    }
    
    np<-length(PFP$pathway_id)
    mutGenes<-rep('',length(np))
    for (i in 1:np){
        currentGenes<-PATH_COLLECTION$HGNC_SYMBOL[[PFP$pathway_id[i]]]
        currentGenes<-intersect(currentGenes,rownames(EM))
        
        if (length(currentGenes)>1){
            freqs<-sort(100*rowSums(sign(EM[currentGenes,]))/ncol(EM),decreasing=TRUE)
            freqs<-freqs[freqs>0]
            mutGenes[i]<-paste(paste(names(freqs),' (',format(freqs,digits=2),' %)',sep=''),collapse=' ')
        }else{
            freqs<-100*sum(sign(EM[currentGenes,]))/ncol(EM)
            mutGenes[i]<-paste(paste(currentGenes,' (',format(freqs,digits=2),' %)',sep=''),collapse=' ')
        }
    }
    
    totres<-cbind(NAMES,
                  PFP$pathway_id,
                  PATH_COLLECTION$Ngenes[PFP$pathway_id],
                  TOTexlength,
                  PFP$pathway_mus,
                  rowSums(PFP$pathway_EM),
                  PFP$pathway_logOddRatios,
                  PFP$pathway_pvals,
                  PFP$pathway_perc_fdr,
                  PFP$pathway_exclusiveCoverage,
                  mutGenes,
                  PFP$pathway_Probability)
    
    colnames(totres)<-c('Pathway',
                        'Internal Id',
                        'n.genes.in.rep',
                        'total Exonic Content Block length',
                        '(E)xpected n.mut in path',
                        '(O)bserved n.mut in path',
                        'log10 (oddRatio=(O/E))',
                        'pval',
                        '% FDR',
                        'ExclusiveCoverage',
                        'Mutated Genes',paste('prob ',colnames(PFP$pathway_Probability)))
    
    write.csv(totres,file=filename,quote=FALSE,row.names=FALSE)
}
SLAPE.heuristic_mut_ex_sorting<-function(mutPatterns){
    
    mutPatterns<-sign(mutPatterns)
    
    if(dim(mutPatterns)[1]==1){
        mutPatterns<-matrix(c(mutPatterns[,order(mutPatterns,decreasing=TRUE)]),
                            1,ncol(mutPatterns),
                            dimnames = list(rownames(mutPatterns),colnames(mutPatterns)))
        
        return(mutPatterns)
    }
    
    if(dim(mutPatterns)[2]==1){
        mutPatterns<-matrix(c(mutPatterns[order(mutPatterns,decreasing=TRUE),]),
                            nrow(mutPatterns),1,
                            dimnames = list(rownames(mutPatters),colnames(mutPatterns)))
        return(mutPatterns)
    }
    
    nsamples<-ncol(mutPatterns)
    
    coveredGenes<-NA
    uncoveredGenes<-rownames(mutPatterns)
    
    if (length(uncoveredGenes)>1){
        
        idNull<-which(colSums(mutPatterns)==0)
        nullCol<-matrix(c(mutPatterns[,idNull]),nrow(mutPatterns),length(idNull),dimnames = list(rownames(mutPatterns),colnames(mutPatterns)[idNull]))
        
        idNonNull<-which(colSums(mutPatterns)>0)
        mutPatterns<-matrix(c(mutPatterns[,idNonNull]),nrow(mutPatterns),length(idNonNull),dimnames=list(rownames(mutPatterns),colnames(mutPatterns)[idNonNull]))
        
        coveredSamples<-NA
        uncoveredSamples<-colnames(mutPatterns)
        BS<-NA
        
        while(length(uncoveredGenes)>0 & length(uncoveredSamples)>0){
            
            patterns<-matrix(c(mutPatterns[uncoveredGenes,uncoveredSamples]),
                             nrow = length(uncoveredGenes),
                             ncol = length(uncoveredSamples),
                             dimnames = list(uncoveredGenes,uncoveredSamples))
            
            if(length(uncoveredGenes)>1){
                bestInClass<-SLE.findBestInClass(patterns)
            }else{
                bestInClass<-uncoveredGenes
            }
            
            if(is.na(BS[1])){
                BS<-bestInClass
            }else{
                BS<-c(BS,bestInClass)
            }
            
            if(is.na(coveredGenes[1])){
                coveredGenes<-bestInClass
            }else{
                coveredGenes<-c(coveredGenes,bestInClass)
            }
            
            uncoveredGenes<-setdiff(uncoveredGenes,coveredGenes)
            toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),nrow = 1,ncol=ncol(patterns),dimnames = list(bestInClass,uncoveredSamples))
            
            if (length(coveredGenes)==1){
                coveredSamples<-names(which(colSums(toCheck)>0))
            }else{
                coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
            }
            
            uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)
            
        }
        
        BS<-c(BS,uncoveredGenes)
        
        CID<-SLE.rearrangeMatrix(mutPatterns,BS)
        
        FINALMAT<-mutPatterns[BS,CID]
        
        FINALMAT<-cbind(FINALMAT,nullCol[rownames(FINALMAT),])
        
        return(FINALMAT)
    }
}
SLAPE.pathvis<-function(EM,PFP,Id,prefName=NULL,PATH='./',PATH_COLLECTION){
    
    genes<-PATH_COLLECTION$HGNC_SYMBOL[[Id]]
    
    nGenesInPath<-length(genes)
    genes<-intersect(genes,rownames(EM))
    
    fn<-paste(PATH,prefName,'_',Id,'_a.pdf',sep='')
    
    toPlot<-matrix(c(EM[genes,]),nrow = length(genes),ncol = ncol(EM),dimnames = list(genes,colnames(EM)))
    
    Pid<-which(PFP$pathway_id==Id)
    
    toPlot<-
        SLAPE.heuristic_mut_ex_sorting(toPlot)
    
    FDR<-PFP$pathway_perc_fdr[which(PFP$pathway_id==Id)]
    
    if (FDR<0){
        FDR<-format(FDR,scientific=TRUE,digits = 3)
    }else{
        FDR<-format(FDR,digits = 3)
    }
    
    NAME<-PATH_COLLECTION$PATHWAY[Id]
    
    NAME<-paste(str_trim(unlist(str_split(NAME,'//'))),collapse='\n')
    
    MAIN<-paste(NAME,'\n\n','SLAPenrich FDR ',FDR,' %',sep='')
    
    pheatmap(toPlot,cluster_rows = FALSE,cluster_cols = FALSE,
             filename = fn,col=c('white','blue'),
             cellheight = 40,
             legend_breaks=c(0,1),
             legend_labels=c('wt','mut'),
             main=MAIN)
    
    pdf(paste(PATH,prefName,'_',PFP$pathway_id[Pid],'_bars.pdf',sep=''))
    
    par(mfrow=c(3,1))
    
    barplot(colSums(EM),las=2,main='n. mutated genes per sample',names.arg = '')
    
    
    barplot(PFP$pathway_Probability[Pid,],las=2,main='p(one mutation in this pathway)',log = 'y',
            xlab=paste('Expected number of samples with mutation in this pathway =',
                       format(sum(PFP$pathway_Probability[Pid,]),digits=2)),names.arg = '')
    
    if(length(genes)>1){
        EM<-as.numeric(sign(colSums(EM[genes,])))
    }else{
        EM<-as.numeric(EM[genes,])
    }
    
    barplot(EM,las=2,
            main=paste('Samples with mutations in this pathway =',sum(EM)),names.arg = '',xlab='samples',yaxt='n')
    
    dev.off()
}
SLAPE.serialPathVis<-function(EM,PFP,fdrth=5,exCovTh=50,PATH='./',PATH_COLLECTION){
    Ids<-which(PFP$pathway_perc_fdr<fdrth & PFP$pathway_exclusiveCoverage>exCovTh)
    
    if (length(Ids)==0){
        warning(paste('No significant enrichments identified at the selected fdr threshold! min fdr % =',
                      min(PFP$pathway_perc_fdr)))
    }else{
        
        print('Producing plots for significantly SL enriched pathways...')
        pb <- txtProgressBar(min=1,max=length(Ids),style=3)
        
        for (i in 1:length(Ids)){
            
            setTxtProgressBar(pb, i)
            SLAPE.pathvis(EM = EM,PFP = PFP,prefName = i,Id = PFP$pathway_id[Ids[i]],PATH = PATH,PATH_COLLECTION = PATH_COLLECTION)
        }
        Sys.sleep(1)
        close(pb)
        print('+ Done!')
    }
}
SLAPE.core_components<-function(PFP,EM,PATH='./',fdrth=Inf,exclcovth=0,PATH_COLLECTION){
    
    filename<-PATH
    EM<-sign(EM)
    
    id<-which(PFP$pathway_exclusiveCoverage>exclcovth & PFP$pathway_perc_fdr<fdrth)
    
    PFP$pathway_id<-PFP$pathway_id[id]
    PFP$pathway_mus<-PFP$pathway_mus[id]
    PFP$pathway_EM<-PFP$pathway_EM[id,]
    PFP$pathway_logOddRatios<-PFP$pathway_logOddRatios[id]
    PFP$pathway_pvals<-PFP$pathway_pvals[id]
    PFP$pathway_perc_fdr<-PFP$pathway_perc_fdr[id]
    PFP$pathway_exclusiveCoverage<-PFP$pathway_exclusiveCoverage[id]
    PFP$pathway_Probability<-PFP$pathway_Probability[id,]
    
    NAMES<-PATH_COLLECTION$PATHWAY[PFP$pathway_id]
    
    NAMES<-str_replace_all(NAMES,',','//')
    
    np<-length(PFP$pathway_id)
    
    print('Producing heatmaps for core-components of enriched pathways...')
    
    
    for (i in 1:np){
        currentGenes<-PATH_COLLECTION$HGNC_SYMBOL[[PFP$pathway_id[i]]]
        currentGenes<-intersect(currentGenes,rownames(EM))
        
        if (i == 1){
            mutG<-currentGenes
        }else{
            mutG<-union(mutG,currentGenes)
        }
    }
    
    mutG<-sort(mutG)
    
    FREQS<-sort(rowSums(EM[mutG,]),decreasing=TRUE)
    
    MM<-SLE.buildPathMemb(mutG,PFP$pathway_id,PATH_COLLECTION)
    
    cm<-fastgreedy.community(graph.incidence(MM))
    
    cmcardinality<-sort(summary(as.factor(cm$membership)),decreasing=TRUE)
    
    communities<-as.numeric(names(cmcardinality))
    ncomm<-length(communities)
    pb <- txtProgressBar(min=1,max=ncomm,style=3)
    
    for (i in 1:ncomm){
        setTxtProgressBar(pb, i)
        
        filen<-paste(filename,'cm_',i,'.pdf',sep='')
        elements<-cm$names[which(cm$membership==communities[i])]
        subData<-MM[intersect(elements,rownames(MM)),intersect(elements,colnames(MM))]
        
        if (length(dim(subData))>0){
            subData<-subData[order(rowSums(subData)),]
            subData<-subData[,order(colSums(subData))]
            verdata<-colnames(subData)
        }else{
            verdata<-intersect(elements,colnames(MM))
        }
        
        verdata<-match(verdata,PFP$pathway_id)
        verdata<- PFP$pathway_perc_fdr[verdata]
        if (length(dim(subData))>0){
            names(verdata)<-PATH_COLLECTION$PATHWAY[as.numeric(colnames(subData))]
            
            SLE.plotMyHeat(subData,100*rowSums(EM[rownames(subData),])/ncol(EM),verdata,filename=filen)
        }else{
            names(verdata)<-PATH_COLLECTION$PATHWAY[as.numeric(intersect(elements,colnames(MM)))]
            SLE.plotMyHeat(matrix(subData,length(subData),1,dimnames = list(names(subData),NULL)),100*rowSums(EM[names(subData),])/ncol(EM),verdata,
                           filename=filen)
        }
        
        
        genesInCoreComponent<-100*rowSums(EM[rownames(subData),])/ncol(EM)
        pathwaysAndEnrichsFDR<-verdata
        
        DATAtoPlot<-list(genesInCoreComponents=genesInCoreComponent,pathwaysAndFDR=pathwaysAndEnrichsFDR)
        
        save(DATAtoPlot,file=paste(filen,'.rdata',sep=''))
        
    }
    Sys.sleep(1)
    close(pb)
}

SLAPE.diff_SLAPE_analysis<-function(EM,contrastMatrix,positiveCondition,negativeCondition,
                                    show_progress=TRUE,display=TRUE,
                                    correctionMethod='fdr',path_probability='Bernoulli',
                                    NSAMPLES=1,NGENES=1,accExLength=TRUE,
                                    BACKGROUNDpopulation=NULL,GeneLenghts,
                                    PATH_COLLECTION,SLAPE.FDRth=5,PATH='./'){
    
    
    
    if(length(positiveCondition)==1){
        positiveSamples<-names(which(contrastMatrix[,positiveCondition]==1))    
    }else{
        positiveSamples<-names(which(rowSums(contrastMatrix[,positiveCondition])==1))
        positiveCondition<-paste(positiveCondition,collapse=', ')
    }
    
    if(length(negativeCondition)==1){
        negativeSamples<-names(which(contrastMatrix[,negativeCondition]==1))
        
    }else{
        negativeSamples<-names(which(rowSums(contrastMatrix[,negativeCondition])==1))
        negativeCondition<-paste(negativeCondition,collapse=', ')
    }
    
    
    positiveSamples<-intersect(positiveSamples,colnames(EM))
    negativeSamples<-intersect(negativeSamples,colnames(EM))
    
    positiveEM<-EM[,positiveSamples]
    negativeEM<-EM[,negativeSamples]
    
    print('Analizing positive population...')
    positive_PFP<-
        SLAPE.analyse(EM = positiveEM,path_probability = path_probability,
                      GeneLenghts = GeneLenghts,
                      show_progress = show_progress,
                      NSAMPLES = 0,
                      NGENES = 0,
                      BACKGROUNDpopulation = BACKGROUNDpopulation,
                      PATH_COLLECTION = PATH_COLLECTION,
                      correctionMethod = correctionMethod)
    print('Done')
    
    print('Analizing negative population...')
    negative_PFP<-
        SLAPE.analyse(EM = negativeEM,path_probability = path_probability,
                      show_progress = show_progress,GeneLenghts = GeneLenghts,
                      NSAMPLES = 0,
                      NGENES = 0,
                      BACKGROUNDpopulation = BACKGROUNDpopulation,
                      PATH_COLLECTION = PATH_COLLECTION,
                      correctionMethod = correctionMethod)
    print('Done')
    
    
    positive.enrichedP.id<-positive_PFP$pathway_id[which(positive_PFP$pathway_perc_fdr<SLAPE.FDRth)]
    negative.enrichedP.id<-negative_PFP$pathway_id[which(negative_PFP$pathway_perc_fdr<SLAPE.FDRth)]
    
    allIDS<-union(positive.enrichedP.id,negative.enrichedP.id)
    
    positive.pvals<-positive_PFP$pathway_pvals[match(allIDS,positive_PFP$pathway_id)]
    negative.pvals<-negative_PFP$pathway_pvals[match(allIDS,negative_PFP$pathway_id)]
    
    names(positive.pvals)<-allIDS
    names(negative.pvals)<-allIDS
    
    differentialEnrichScores<- -log10(positive.pvals)+log10(negative.pvals)
    differentialEnrichScores<- sort(differentialEnrichScores,decreasing=TRUE)
    
    ids<-names(differentialEnrichScores)
    positiveFDR<-
        positive_PFP$pathway_perc_fdr[match(ids,positive_PFP$pathway_id)]
    negativeFDR<-
        negative_PFP$pathway_perc_fdr[match(ids,negative_PFP$pathway_id)]
    
    names(positiveFDR)<-ids
    names(negativeFDR)<-ids
    
    positive_pathway_EM<-positive_PFP$pathway_EM[match(ids,positive_PFP$pathway_id),]
    negative_pathway_EM<-negative_PFP$pathway_EM[match(ids,negative_PFP$pathway_id),]
    
    d <- dist(t(positive_pathway_EM),method = 'euclidean') # distance matrix
    fit <- hclust(d,method = 'complete')
    positive_pathway_EM<-
        positive_pathway_EM[,fit$labels[fit$order]]
    
    d <- dist(t(negative_pathway_EM),method = 'euclidean') # distance matrix
    fit <- hclust(d,method = 'complete')
    negative_pathway_EM<-
        negative_pathway_EM[,fit$labels[fit$order]]
    
    
    COMPD<-cbind(positive_pathway_EM,
                 negative_pathway_EM)
    
    
    FDRs<-cbind(positiveFDR,negativeFDR)
    FDRs<-FDRs[rownames(COMPD),]/100
    
    FDRs<- -log10(FDRs+.Machine$double.eps)
    
    currentNames<-PATHCOM_HUMAN$PATHWAY[as.numeric(rownames(FDRs))]
    suffixes<-rep('',length(currentNames))
    suffixes[which(str_length(currentNames)>50)]<-' ...'
    
    currentNames<-str_sub(currentNames,start = 1,end = 50)
    currentNames<-paste(currentNames,suffixes,sep='')
    rownames(FDRs)<-currentNames
    
    colnames(FDRs)<-c(positiveCondition,negativeCondition)
    
    rownames(COMPD)<-rownames(FDRs)
    
    
    annotation_col = data.frame(SampleType = factor(c(rep(positiveCondition, ncol(positive_pathway_EM)),rep(negativeCondition, ncol(negative_pathway_EM)))))
    rownames(annotation_col)<-colnames(COMPD)
    
    COMPD<-COMPD[c(1:30,(nrow(COMPD)-29):nrow(COMPD)),]
    
    if(display){
        pdf(paste(PATH,'diffPathAlter.pdf',sep=''),10,10)
        pheatmap(COMPD,col=c('white','blue'),annotation_col = annotation_col,show_colnames = FALSE,
                 cluster_rows = FALSE,cluster_cols = FALSE)
        dev.off()
    }
    
    annotation_col = data.frame(CellType = factor(c(positiveCondition,negativeCondition)))
    rownames(annotation_col)<-colnames(FDRs)
    
    
    if(display){
        pdf(paste(PATH,'diffPathAlterHM.pdf',sep=''),10,10)
        pheatmap(FDRs,cluster_rows = FALSE,cluster_cols = FALSE,col=colorRampPalette(colors = c('black','purple'))(100),annotation_col=annotation_col,show_colnames = FALSE)
        dev.off()
        pdf(paste(PATH,'diffPathAlterBP.pdf',sep=''),10,10)
        par(mar=c(4,25,4,4))
        barplot(rev(FDRs[,1]-FDRs[,2]),horiz = TRUE,las=2,xlim=c(-7,7),cex.names = 0.7)
        dev.off()
    }
    
    diffSLenrich<-FDRs[,1]-FDRs[,2]
    
    RES<-cbind(FDRs,diffSLenrich)
    colnames(RES)[3]<-'Diff.SL.Enrich'
    
    return(RES)
}



#not documented

SLAPE.vignette<-function(){
    LungTSPvariants<-read.table('../../Data/SLAPenrich/externalData/supplementary_table_2.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)
    
    GeneIds<-unique(LungTSPvariants$HUGO.Symbol)
    GeneIds<-setdiff(GeneIds,'')
    
    SampleIds<-unique(LungTSPvariants$Tumor.ID)
    SampleIds<-as.character(setdiff(SampleIds,NA))
    
    Dataset<-matrix(0,length(GeneIds),length(SampleIds),dimnames = list(GeneIds,SampleIds))
    
    
    nvariants<-nrow(LungTSPvariants)
    for (i in 1:nvariants){
        sampleId<-as.character(LungTSPvariants$Tumor.ID[i])
        geneId<-LungTSPvariants$HUGO.Symbol[i]
        if (geneId!='' & sampleId!='NA'){
            Dataset[geneId,sampleId]<-Dataset[geneId,sampleId]+1    
        }
    }
    
    ####to change
    load('data/SLAPE.hgnc.table_20160210.Rdata')
    load('data/SLAPE.20160211_MSigDB_KEGG_hugoUpdated.Rdata')
    load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")
    
    
    Dataset<-
        SLAPE.check_and_fix_gs_Dataset(Dataset,updated.hgnc.table = updated.hgnc.table)
    
    PFPw<-
        SLAPE.analyse(EM = Dataset,PATH_COLLECTION = KEGG_PATH,
                      show_progress = TRUE,
                      NSAMPLES = 1,
                      NGENES = 1,
                      accExLength = TRUE,
                      BACKGROUNDpopulation = rownames(Dataset),
                      path_probability = 'Bernoulli',
                      GeneLenghts = GECOBLenghts)
    
    SLAPE.pathvis(EM = Dataset,PFP = PFPw,Id = PFPw$pathway_id[which(PFPw$pathway_exclusiveCoverage>80)[1]],
                  PATH = './',PATH_COLLECTION = KEGG_PATH)
    
    
    SLAPE.write.table(PFP = PFPw,EM = Dataset,filename = "../../RESULTS/SLAPenrich/LungDS_KEGG_enrichments.csv",
                      fdrth=5,exclcovth = 50,PATH_COLLECTION = KEGG_PATH,GeneLenghts = GECOBLenghts)
    
    load('data/caseStudy_clinicalInfos.Rdata')
    load('data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.Rdata')
    
    RES<-
        SLAPE.diff_SLAPE_analysis(EM = Dataset,contrastMatrix = CLINIC_MAT,
                                  BACKGROUNDpopulation = rownames(Dataset),
                                  SLAPE.FDRth = 5,display = TRUE,
                                  positiveCondition = 'SS_CurrentSmoker',
                                  negativeCondition = 'SS_Never',
                                  PATH_COLLECTION = PATHCOM_HUMAN,
                                  PATH = 'temp.results/')
    
    RES1<-SLAPE.diff_SLAPE_analysis(EM = Dataset,contrastMatrix = CLINIC_MAT,
                                    BACKGROUNDpopulation = rownames(Dataset),
                                    SLAPE.FDRth = 5,display = FALSE,
                                    positiveCondition = "BAC_Type_nonMucinous",
                                    negativeCondition = "BAC_Type_Mucinous",
                                    PATH_COLLECTION = PATHCOM_HUMAN,
                                    PATH = 'temp.results/')
    
}

