library(HGNChelper)
library(poibin)
library(stringr)
library(reshape)
library(png)
library(plotrix)
library(ecolitk)

source('R/SLAPenrich.R')
load('data/SLAPE.hgnc.table_20160210.RData')
load('data/SLAPE.20160526_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.rdata')
load('data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160510.RData')
load('data/SLAPE.hgnc.table_20160210.RData')

nGenesPerHallmark<-function(RES,PATHCOM_HUMAN){
    
    uhm<-
        unique(RES$HallMark_Table[,2])
    
    
    for (i in 1:length(uhm)){
        
        currentG<-unique(unlist(PATHCOM_HUMAN$HGNC_SYMBOL[RES$HallMark_Table[which(RES$HallMark_Table[,2]==uhm[i]),1]]))
        names(currentG)<-rep(uhm[i],length(currentG))
        
        if (i == 1){
            Genes<-currentG
        }else{
            Genes<-c(Genes,currentG)   
        }
    }
    
    return(Genes)
}

write.csv.file.hallmarks<-function(RES,path='./'){
    
    npath<-nrow(RES$HallMark_Table)
    
    for (i in 1:npath){
 
        currentIdx<-RES$HallMark_Table[i,1]
        currentHM<-RES$HallMark_Table[i,2]
        currentPT<-as.character(RES$HallMark_Table[,3])[i]
        
        
        id<-match(currentPT,RES$PATHCOLL$PATHWAY)
        
        currentPT<-str_replace_all(currentPT,pattern = '// ',replacement = ' // ')
        currentsource<-RES$PATHCOLL$SOURCE[[id]]
        currentsource<-unlist(str_split(currentsource,'datasource: '))
        currentsource<-currentsource[currentsource!='']
        currentsource<-paste(currentsource, collapse=' // ')
        
        currentGenes<-paste(sort(RES$PATHCOLL$HGNC_SYMBOL[[id]]),collapse=', ')
        
        currentLINE<-c(currentIdx,
                       currentHM,
                       currentPT,
                       currentsource,
                       currentGenes)
        
        if (i == 1){
            restable<-currentLINE
        }else{
            restable<-rbind(restable,currentLINE)
        }
    }
           
    rownames(restable)<-NULL
    colnames(restable)<-c('internal id','Hallmark','Pathway(s)','Source(s)','Genes')
    
    save(restable,file=paste(path,'HallMarks_pathways_Genes.rdata',sep=''))
    write.table(restable,quote=FALSE,sep='\t',row.names = FALSE,file=paste(path,'HallMarks_pathways_Genes.tsv',sep=''))
}


RES<-SLAPE.integrateHallmarks(PATHCOLL = PATHCOM_HUMAN)

write.csv.file.hallmarks(RES,path = 'data/20160719_')


HM_TABLE<-
    RES$HallMark_Table

PATH_COLLECTION<-
    RES$PATHCOLL

DATA_VEC<-
    list(DATA=summary(as.factor(HM_TABLE$Hallmark)),HM=names(summary(as.factor(HM_TABLE$Hallmark))))

GENES<-
    nGenesPerHallmark(RES,PATHCOM_HUMAN)

GENE_VEC<-
    list(DATA=summary(as.factor(names(GENES))),HM=names(summary(as.factor(names(GENES)))))

SLE.radialHallmarkPlot(DATA_VEC,mc = 50,RADIUS = 1,LABEL = 'n. pathways',SCALE = 0.5)
SLE.radialHallmarkPlot(GENE_VEC,mc = 50,RADIUS = 1.7,LABEL = 'n. genes',SCALE = 0.5,TOADD = TRUE)
