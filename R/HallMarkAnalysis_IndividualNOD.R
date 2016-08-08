library(HGNChelper)
library(poibin)
library(stringr)
library(reshape)
library(png)
library(ecolitk)
source('R/SLAPenrich.R')
source('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_CODE/MISC/my.hypTest.R')

load('data/hm_colors.rdata')

TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','PRAD')

TPR<-rep(0,length(TYPES))
PPV<-rep(0,length(TYPES))
COR<-rep(0,length(TYPES))

names(TPR)<-TYPES
names(PPV)<-TYPES
names(COR)<-TYPES

for (TYPE in TYPES){

    
    
#    pdf(paste('../../RESULTS/SLAPenrich/NOD_radialPlots/',TYPE,'.pdf',sep=''),width = 18,height = 5)
#RES<-SLE.HallmarkAnalysis(TYPE=TYPE,PATH_COLLECTION = PATH_COLLECTION,HM_TABLE = RES$HallMark_Table,PATH = '../../RESULTS/SLAPenrich/PT_HM_20160718/',
#                          removeDrivers=TRUE)
    par(mfrow=c(1,4))
    par(xpd=TRUE)
    
    load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPE,'_HM.rdata',sep=''))
    tmp<- (RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
    names(tmp)<-RES$Hallmark
    PATT_VEC_ALL<-list(DATA=tmp,HM=rownames(RES))
    SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_ALL,mc = 50, RADIUS = 3,SCALE =1.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE,
                           title='pathways')
    
    
    text(-7,6,TYPE,cex=4,pos = 4)
    
    pt<-unique(RES$Idx)
    id<-match(pt,RES$Idx)
    
    N<-length(id)
    RES<-RES[id,]
    
    enriched_all<-which(RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
    n<-length(enriched_all)
    
    ttmp<-tmp
    un<-unique(names(tmp))
    for (i in 1:length(un)){
        ttmp[which(names(tmp)==un[i])]<-rep(sum(tmp[which(names(tmp)==un[i])])/length(which(names(tmp)==un[i])),
                                            length(which(names(tmp)==un[i])))
    }
    PATT_VEC_ALL_HM<-list(DATA=ttmp,HM=rownames(RES))
        
    
    load(paste('../../RESULTS/SLAPenrich/PT_HM_20160718/',TYPE,'_HM.rdata',sep=''))
    tmp<- (RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
    names(tmp)<-RES$Hallmark
    PATT_VEC_NOD<-list(DATA=tmp,HM=rownames(RES))
    SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_NOD,mc = 50, RADIUS =2,SCALE = 1,hanna = FALSE,scaleNumb=FALSE,resort = FALSE,
                           TOADD = TRUE,
                           title='')
    
    SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_ALL_HM,maxscale=1,mc = 50,title = 'hallmarks',
                           RADIUS = 3,TOADD = FALSE,SCALE = 1.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)
    
    
    ttmp<-tmp
    un<-unique(names(tmp))
    for (i in 1:length(un)){
        ttmp[which(names(tmp)==un[i])]<-rep(sum(tmp[which(names(tmp)==un[i])])/length(which(names(tmp)==un[i])),
                                            length(which(names(tmp)==un[i])))
    }
    PATT_VEC_NOD_HM<-list(DATA=ttmp*5,HM=rownames(RES))
    SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_NOD_HM,maxscale=2,mc = 50, RADIUS = 2,
                           TOADD =TRUE,SCALE = 1,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)
    
    
    
    pt<-unique(RES$Idx)
    id<-match(pt,RES$Idx)
    
    RES<-RES[id,]
    
    enriched_nod<-which(RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
    k<-length(enriched_nod)
    
    x<-length(intersect(enriched_all,enriched_nod))
    
    print(sum(PATT_VEC_ALL$DATA*PATT_VEC_NOD$DATA)/sum(PATT_VEC_ALL$DATA)*100)
    print(sum(PATT_VEC_ALL$DATA*PATT_VEC_NOD$DATA)/sum(PATT_VEC_NOD$DATA)*100)
    
    p<-my.hypTest(x,k,n,N)
    
    ALL_HM<-PATT_VEC_ALL_HM$DATA[unique(names(PATT_VEC_ALL_HM$DATA))]
    NOD_HM<-PATT_VEC_NOD_HM$DATA[unique(names(PATT_VEC_NOD_HM$DATA))]/5
    
    
    TPR[TYPE]<-100*x/n
    PPV[TYPE]<-100*x/k
    
    barplot(100*c(x/n,x/k),col='gray',las=2,srt = 45,ylim=c(0,100),
            ylab='% pathways',main='novel vs. all variant analysis',cex.axis = 1.5,
            cex.lab = 1.5,cex.main=1.5)
    text(c(0.75,1.95), par("usr")[3]-5, 
         srt = 60, adj= 1, xpd = TRUE,
         labels = c('TPR','PPV'), cex=2)
    
    legend('topleft',legend=paste('p = ',format(p,scientific=TRUE,digits=2),sep=''),cex = 1.5,bty = 'n')
    
    COR[TYPE]<-cor(ALL_HM,NOD_HM,method = 'spearman')
    
    plot(ALL_HM,NOD_HM,col=hm_colors[names(ALL_HM)],xlim=c(0,max(ALL_HM)),ylim=c(0,max(NOD_HM)),pch=16,cex=2,
         cex.lab=1.5,cex.axis=1.5,xlab='All variants',ylab='Novel variants',cex.main=1.5,
         main=paste('Hallmark footprints R = ',format(cor(ALL_HM,NOD_HM,method = 'spearman'),digits=2)))
    
 #   dev.off()
}
