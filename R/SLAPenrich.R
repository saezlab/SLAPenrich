# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

SLE.hypTest<-function(x,k,n,N){
    PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
    return(PVALS)
}
SLE.HeuristicMutExSorting<-function(mutPatterns){

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
SLE.buildPathMemb<-function(GENES,pathway_ids){
    np<-length(pathway_ids)
    ngenes<-length(GENES)
    MM<-matrix(0,ngenes,np,dimnames = list(GENES,pathway_ids))

    for (i in 1:np){
        MM[,i]<-is.element(GENES,PATHCOM_HUMAN$HGNC_SYMBOL[[pathway_ids[i]]])+0
    }

    return(MM)
}
SLE.KnownGenesInRefGenome<-function(){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    gn <- sort(genes(txdb))
    gn <- unlist(gn$gene_id)
    return(gn)
}
SLE.computeGeneLengths<-function(){

    exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')

    egs    = KnownGenesInRefGenome()
    egs<-intersect(egs,keys(org.Hs.egACCNUM))

    GL<-rep(NA,length(egs))
    names(GL)<-egs

    GL<-sapply(egs,function(eg)
    {
        print(which(egs==eg))

        exons = exons.db[[eg]]
        if(length(exons)>0){
            exons = reduce(exons)
            L<-sum( width(exons) )
        }else{
            L<-NA
        }
    })

    names(GL)<-getSYMBOL(names(GL), data='org.Hs.eg')
    return(GL)
}
SLE.assignTotalLengthTOpathways<-function(){
    np<-length(PATHCOM_HUMAN$PATHWAY)

    Glenghts<-list()
    TOTAL_G_LENGHT<-vector()
    for (i in 1:np){
        print(i)
        Glenghts[[i]]<-gene_in_path_lengths[PATHCOM_HUMAN$HGNC_SYMBOL[[i]]]
        TOTAL_G_LENGHT[i]<-sum(gene_in_path_lengths[PATHCOM_HUMAN$HGNC_SYMBOL[[i]]],na.rm = TRUE)
    }

    names(Glenghts)<-PATHCOM_HUMAN$PATHWAY
    names(TOTAL_G_LENGHT)<-PATHCOM_HUMAN$PATHWAY

    PATHCOM_HUMAN$Glengths<-Glenghts
    PATHCOM_HUMAN$TOTAL_G_LENGHT<-TOTAL_G_LENGHT

    return(PATHCOM_HUMAN)

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


SLAPE.readDataset<-function(filename){
    fc<-as.matrix(read.csv(filename,row.names=1))
    fc[is.na(fc)]<-0

    return(fc)
}
SLAPE.Analyse<-function(wBEM,show_progress=TRUE,correctionMethod='fdr',NSAMPLES=1,NGENES=1,accExLength=TRUE,
                                     BACKGROUNDpopulation=NULL){

    if(length(BACKGROUNDpopulation)>0){
            if(length(duplicated(BACKGROUNDpopulation))>0){
                warning('Duplicated gene names found in the background population')
                }
            BACKGROUNDpopulation<-unique(BACKGROUNDpopulation)
        }

    if(sum(is.na(wBEM))>0){
        error('The inputted BEM contains NA!')
    }

    if(accExLength){
        cat('Sample fingerprinting for pathawy alterations (accounting for gene exonic length)...\n')
        gLenghts<-unlist(PATHCOM_HUMAN$Glengths)
        if(!length(BACKGROUNDpopulation)){
            N<-sum(gLenghts[unique(names(gLenghts))],na.rm = TRUE)
        }else{
            N<-all_genes_exonic_lengths[unique(BACKGROUNDpopulation)]
        }
    }else{
        cat('Sample fingerprinting for pathawy alterations...\n')
        if(!length(BACKGROUNDpopulation)){
            N<-length(unique(unlist(PATHCOM_HUMAN$Glengths)))
        }else{
            N<-length(BACKGROUNDpopulation)
            }
        }

    np<-length(PATHCOM_HUMAN$PATHWAY)

    nsamples<-ncol(wBEM)

    PN<-PATHCOM_HUMAN$PATHWAY

    pathway_BEM<-matrix(0,nrow=np,ncol=nsamples,dimnames=list(1:np,colnames(wBEM)))
    pathway_Probability<-matrix(0,nrow=np,ncol=nsamples,dimnames=list(1:np,colnames(wBEM)))
    pathway_grandTotal<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_pvals<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_mus<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathwayExclusive_coverage<-matrix(0,nrow=1,ncol=np,dimnames=list(NULL,1:np))
    pathway_individualBEMs<-list()

    if(show_progress){
        pb <- txtProgressBar(min=1,max=nsamples*np,style=3)
    }

    toExclude<-rep(FALSE,np)
    for (i in 1:np){

        currentGeneSet<-intersect(PATHCOM_HUMAN$HGNC_SYMBOL[[i]],rownames(wBEM))

        GLENGHTS<-PATHCOM_HUMAN$Glengths[[i]][currentGeneSet]

        if (length(which(!is.na(GLENGHTS)))==0 | length(currentGeneSet)<=NGENES){
            toExclude[i]<-TRUE
        }

        if(accExLength){
            n <- sum(PATHCOM_HUMAN$Glengths[[i]][currentGeneSet],na.rm = TRUE)
        }else{
            n <-length(currentGeneSet)
        }


        if (length(currentGeneSet)>1){
            pathway_individualBEMs[[i]]<-wBEM[currentGeneSet,]
        }else{
            pathway_individualBEMs[[i]]<-matrix(wBEM[currentGeneSet,],1,ncol(wBEM),dimnames = list(currentGeneSet,colnames(wBEM)))
        }

        sampleCovered<-sum(colSums(pathway_individualBEMs[[i]])>0)
        exclusiveSampleCovered<-sum(colSums(pathway_individualBEMs[[i]])==1)

        pathwayExclusive_coverage[i]<-100*exclusiveSampleCovered/sampleCovered

        if (length(currentGeneSet)>0){
            pathway_grandTotal[i]<-sum(unlist(c(wBEM[currentGeneSet,])))
        }else{
            pathway_grandTotal[i]<-0
        }

        for (j in 1:nsamples){
            k <- sum(wBEM[,j])
            x<-1
            if(length(currentGeneSet)>0){
                pathway_BEM[i,j]<-sign(sum(wBEM[currentGeneSet,j]))
            }
            pathway_Probability[i,j]<-SLE.hypTest(x,k,n,N)
            if(show_progress){setTxtProgressBar(pb, (i-1)*nsamples+j)}
        }
    }

    if(show_progress){
        Sys.sleep(1)
        close(pb)
    }

    pathway_id<-1:length(PATHCOM_HUMAN$PATHWAY)

    pathway_mus<-rowSums(pathway_Probability)
    pathway_var<-rowSums(1-pathway_Probability)*pathway_Probability

    pathway_logOddRatios<-log10(rowSums(pathway_BEM)/pathway_mus)

    for (i in 1:nrow(pathway_Probability)){
        pathway_pvals[i]<-sum(dpoibin(rowSums(pathway_BEM)[i]:ncol(wBEM),pathway_Probability[i,]))
    }

    id<-which(rowSums(pathway_BEM)>=NSAMPLES)

    pathway_id<-pathway_id[id]
    pathway_BEM<-pathway_BEM[id,]
    pathway_Probability<-pathway_Probability[id,]
    pathway_pvals<-pathway_pvals[id]
    pathway_mus<-pathway_mus[id]
    pathway_logOddRatios<-pathway_logOddRatios[id]
    pathway_individualBEMs<-pathway_individualBEMs[id]
    pathwayExclusive_coverage<-pathwayExclusive_coverage[id]
    toExclude<-toExclude[id]

    pathway_id<-pathway_id[order(pathway_pvals)]
    pathway_BEM<-pathway_BEM[order(pathway_pvals),]
    pathway_Probability<-pathway_Probability[order(pathway_pvals),]
    pathway_logOddRatios<-pathway_logOddRatios[order(pathway_pvals)]
    pathway_mus<-pathway_mus[order(pathway_pvals)]
    pathway_individualBEMs<-pathway_individualBEMs[order(pathway_pvals)]
    pathwayExclusive_coverage<-pathwayExclusive_coverage[order(pathway_pvals)]
    toExclude<-toExclude[order(pathway_pvals)]
    pathway_pvals<-sort(pathway_pvals)

    pathway_id<-pathway_id[which(!toExclude)]
    pathway_BEM<-pathway_BEM[which(!toExclude),]
    pathway_Probability<-pathway_Probability[which(!toExclude),]
    pathway_logOddRatios<-pathway_logOddRatios[which(!toExclude)]
    pathway_mus<-pathway_mus[which(!toExclude)]
    pathway_individualBEMs<-pathway_individualBEMs[which(!toExclude)]
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

    return(list(pathway_id=pathway_id,
                pathway_BEM=pathway_BEM,
                pathway_Probability=pathway_Probability,
                pathway_mus=pathway_mus,
                pathway_logOddRatios=pathway_logOddRatios,
                pathway_pvals=pathway_pvals,
                pathway_perc_fdr=pathway_perc_fdr,
                pathway_exclusiveCoverage=pathwayExclusive_coverage,
                pathway_individualBEMs=pathway_individualBEMs))
}
SLAPE.write.table<-function(PFP,BEM,filename='',fdrth=Inf,exclcovth=0){

        id<-which(PFP$pathway_exclusiveCoverage>exclcovth & PFP$pathway_perc_fdr<fdrth)

        PFP$pathway_id<-PFP$pathway_id[id]
        PFP$pathway_mus<-PFP$pathway_mus[id]
        PFP$pathway_BEM<-PFP$pathway_BEM[id,]
        PFP$pathway_logOddRatios<-PFP$pathway_logOddRatios[id]
        PFP$pathway_pvals<-PFP$pathway_pvals[id]
        PFP$pathway_perc_fdr<-PFP$pathway_perc_fdr[id]
        PFP$pathway_exclusiveCoverage<-PFP$pathway_exclusiveCoverage[id]
        PFP$pathway_Probability<-PFP$pathway_Probability[id,]

        NAMES<-PATHCOM_HUMAN$PATHWAY[PFP$pathway_id]

        NAMES<-str_replace_all(NAMES,',','//')

        TOTexlength<-unlist(lapply(PATHCOM_HUMAN$Glengths[PFP$pathway_id],FUN = 'sum',na.rm=TRUE))

        np<-length(PFP$pathway_id)
        mutGenes<-rep('',length(np))
        for (i in 1:np){
            currentGenes<-PATHCOM_HUMAN$HGNC_SYMBOL[[PFP$pathway_id[i]]]
            currentGenes<-intersect(currentGenes,rownames(BEM))

            if (length(currentGenes)>1){
                freqs<-sort(100*rowSums(sign(BEM[currentGenes,]))/ncol(BEM),decreasing=TRUE)
                freqs<-freqs[freqs>0]
                mutGenes[i]<-paste(paste(names(freqs),' (',format(freqs,digits=2),' %)',sep=''),collapse=' ')
            }else{
                freqs<-100*sum(sign(BEM[currentGenes,]))/ncol(BEM)
                mutGenes[i]<-paste(paste(currentGenes,' (',format(freqs,digits=2),' %)',sep=''),collapse=' ')
            }
        }

        totres<-cbind(NAMES,
                      PFP$pathway_id,
                      PATHCOM_HUMAN$Ngenes[PFP$pathway_id],
                      TOTexlength,
                      PFP$pathway_mus,
                      rowSums(PFP$pathway_BEM),
                      PFP$pathway_logOddRatios,
                      PFP$pathway_pvals,
                      PFP$pathway_perc_fdr,
                      PFP$pathway_exclusiveCoverage,
                      mutGenes,
                      PFP$pathway_Probability)

        colnames(totres)<-c('Pathway','Internal Id','n.genes.in.rep',
                            'total Exonic length',
                            '(E)xpected n.mut in path',
                            '(O)bserved n.mut in path','log10 (oddRatio=(O/E))',
                            'pval','% FDR','ExclusiveCoverage','Mutated Genes',paste('prob ',colnames(PFP$pathway_Probability)))

        write.csv(totres,file=filename,quote=FALSE,row.names=FALSE)
    }
SLAPE.pathVis<-function(BEM,PFP,Id,i=NULL,PATH='./'){

    genes<-PATHCOM_HUMAN$HGNC_SYMBOL[[Id]]

    nGenesInPath<-length(genes)
    genes<-intersect(genes,rownames(BEM))

    fn<-paste(PATH,i,'_',Id,'_a.pdf',sep='')

    toPlot<-matrix(c(BEM[genes,]),nrow = length(genes),ncol = ncol(BEM),dimnames = list(genes,colnames(BEM)))

    Pid<-which(PFP$pathway_id==Id)

    toPlot<-SLE.HeuristicMutExSorting(toPlot)

    FDR<-PFP$pathway_perc_fdr[which(PFP$pathway_id==Id)]

    if (FDR<0){
        FDR<-format(FDR,scientific=TRUE,digits = 3)
    }else{
        FDR<-format(FDR,digits = 3)
    }

    NAME<-PATHCOM_HUMAN$PATHWAY[Id]

    NAME<-paste(str_trim(unlist(str_split(NAME,'//'))),collapse='\n')

    MAIN<-paste(NAME,'\n\n','SLAPenrich FDR ',FDR,' %',sep='')

    pheatmap(toPlot,cluster_rows = FALSE,cluster_cols = FALSE,
             filename = fn,col=c('white','blue'),
             cellheight = 40,
             legend_breaks=c(0,1),
             legend_labels=c('wt','mut'),
             main=MAIN)

    pdf(paste(PATH,i,'_',PFP$pathway_id[Pid],'_bars.pdf',sep=''))

    par(mfrow=c(3,1))

    barplot(colSums(BEM),las=2,main='n. mutated genes per sample',names.arg = '')


    barplot(PFP$pathway_Probability[Pid,],las=2,main='p(one mutation in this pathway)',log = 'y',
            xlab=paste('Expected number of samples with mutation in this pathway =',
                       format(sum(PFP$pathway_Probability[Pid,]),digits=2)),names.arg = '')

    if(length(genes)>1){
        bem<-as.numeric(sign(colSums(BEM[genes,])))
    }else{
        bem<-as.numeric(BEM[genes,])
    }

    barplot(bem,las=2,
            main=paste('Samples with mutations in this pathway =',sum(bem)),names.arg = '',xlab='samples',yaxt='n')

    dev.off()
}
SLAPE.serialPathVis<-function(BEM,PFP,fdrth=5,exCovTh=50,PATH='./'){
    Ids<-which(PFP$pathway_perc_fdr<fdrth & PFP$pathway_exclusiveCoverage>exCovTh)

    if (length(Ids)==0){
        warning(paste('No significant enrichments identified at the selected fdr threshold! min fdr % =',
                      min(PFP$pathway_perc_fdr)))
    }else{

        print('Producing plots for significantly SL enriched pathways...')
        pb <- txtProgressBar(min=1,max=length(Ids),style=3)

        for (i in 1:length(Ids)){

            setTxtProgressBar(pb, i)
            SLAPE.pathVis(BEM = BEM,PFP = PFP,i = i,Id = PFP$pathway_id[Ids[i]],PATH = PATH)
        }
        Sys.sleep(1)
        close(pb)
        print('+ Done!')
    }
}
SLAPE.coreComponents<-function(PFP,BEM,filename='',fdrth=Inf,exclcovth=0){

    BEM<-sign(BEM)

    id<-which(PFP$pathway_exclusiveCoverage>exclcovth & PFP$pathway_perc_fdr<fdrth)

    PFP$pathway_id<-PFP$pathway_id[id]
    PFP$pathway_mus<-PFP$pathway_mus[id]
    PFP$pathway_BEM<-PFP$pathway_BEM[id,]
    PFP$pathway_logOddRatios<-PFP$pathway_logOddRatios[id]
    PFP$pathway_pvals<-PFP$pathway_pvals[id]
    PFP$pathway_perc_fdr<-PFP$pathway_perc_fdr[id]
    PFP$pathway_exclusiveCoverage<-PFP$pathway_exclusiveCoverage[id]
    PFP$pathway_Probability<-PFP$pathway_Probability[id,]

    NAMES<-PATHCOM_HUMAN$PATHWAY[PFP$pathway_id]

    NAMES<-str_replace_all(NAMES,',','//')

    np<-length(PFP$pathway_id)

    print('Producing heatmaps for core-components of enriched pathways...')


    for (i in 1:np){
        currentGenes<-PATHCOM_HUMAN$HGNC_SYMBOL[[PFP$pathway_id[i]]]
        currentGenes<-intersect(currentGenes,rownames(BEM))

        if (i == 1){
            mutG<-currentGenes
        }else{
            mutG<-union(mutG,currentGenes)
        }
    }

    mutG<-sort(mutG)

    FREQS<-sort(rowSums(BEM[mutG,]),decreasing=TRUE)

    MM<-SLE.buildPathMemb(mutG,PFP$pathway_id)

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
            names(verdata)<-PATHCOM_HUMAN$PATHWAY[as.numeric(colnames(subData))]

            SLE.plotMyHeat(subData,100*rowSums(BEM[rownames(subData),])/ncol(BEM),verdata,filename=filen)
        }else{
            names(verdata)<-PATHCOM_HUMAN$PATHWAY[as.numeric(intersect(elements,colnames(MM)))]
            SLE.plotMyHeat(matrix(subData,length(subData),1,dimnames = list(names(subData),NULL)),100*rowSums(BEM[names(subData),])/ncol(BEM),verdata,
                       filename=filen)
        }

        
        genesInCoreComponent<-100*rowSums(BEM[rownames(subData),])/ncol(BEM)
        pathwaysAndEnrichsFDR<-verdata

        DATAtoPlot<-list(genesInCoreComponents=genesInCoreComponent,pathwaysAndFDR=pathwaysAndEnrichsFDR)

        save(DATAtoPlot,file=paste(filen,'.rdata',sep=''))

    }
    Sys.sleep(1)
    close(pb)
}

