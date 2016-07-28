source('R/HallMarkAnalysis_Preamble.R')
load('data/hm_colors.rdata')

fixDriverSymbols<-function(GL){
    tmp<-checkGeneSymbols(GL)$Suggested.Symbol
    tmp<-tmp[!is.na(tmp)]
    return(tmp)
}

TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV','PRAD')

load('data/IntoGen_cancer_drivers_13_06_2014.rdata')



for (TYPE in TYPES){
    print(TYPE)
    currentDrivers<-
        cancerDrivers$ActingDriver_Symbol[which(cancerDrivers$Tumor_Type==TYPE)]
        
    currentDrivers<-
        fixDriverSymbols(currentDrivers)
    
    load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPE,'_HM.rdata',sep=''))
    BG<-unique(unlist(PATHCOM_HUMAN$HGNC_SYMBOL[RES$Idx]))
    
    
    notMatched<-
        setdiff(currentDrivers,BG)
    
    currentDrivers<-
        intersect(currentDrivers,BG)
        
    
    tmp<- which(RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
    
    RES<-RES[tmp,]
    
    um<-c("Activating Invasion and Metastasis",
          "Avoiding Immune Destruction",
          "Deregulating Cellular Energetics",
          "Enabling Replicative Immortality",
          "Evading Growth Suppressors",
          "Genome Instability and Mutation",
          "Inducing Angiogenesis",
          "Resisting Cell Death",
          "Sustaining Proliferative Signaling",
          "Tumour-Promoting Inflammation")
    
    GM<-matrix(0,length(currentDrivers),length(um),dimnames = list(currentDrivers,um))
 
    for (i in 1:nrow(RES)){
        GS<-
            PATHCOM_HUMAN$HGNC_SYMBOL[[RES$Idx[i]]]
        
        GS<-intersect(GS,currentDrivers)
        
        GM[GS,RES$Hallmark[i]]<-1
    }
       
    nNotMapped<-length(notMatched)
    nMultipleHM<-length(which(rowSums(GM)>1))
    
    idx<-which(rowSums(GM)==1)
    GM<-GM[idx,]
    
    remaining<-colSums(GM)
    notRecovered<-length(currentDrivers)-sum(c(nMultipleHM,remaining))
    
    cr<-matrix(c(remaining,nMultipleHM,notRecovered,nNotMapped),nrow = 13,ncol=1,
           dimnames = list(c(names(remaining),'multiple HM','others','not mapped'),TYPE))
    
    if (TYPE=='BRCA'){
        RESTOT<-cr
    }else{
        RESTOT<-cbind(RESTOT,cr)
    }
}

TOTAL<-colSums(RESTOT)
TOTALmapped<-colSums(RESTOT[1:(nrow(RESTOT)-1),])

RESTOT<-RESTOT[1:(nrow(RESTOT)-1),]
norml<-RESTOT/t(matrix(rep(TOTALmapped,nrow(RESTOT)),ncol(RESTOT),nrow(RESTOT)))

norml<-norml[,order(colSums(norml[1:(nrow(norml)-1),]),decreasing=TRUE)]

tmp<-barplot(100*norml[1:(nrow(norml)-1),],col=c(hm_colors[rownames(RESTOT)[1:10]],'white'),ylim=c(0,100),las=2,ylab='% mapped intOGen CGs')

text(tmp,colSums(100*norml[1:(nrow(norml)-1),])+2,paste(TOTALmapped,'/',TOTAL,sep=''))



