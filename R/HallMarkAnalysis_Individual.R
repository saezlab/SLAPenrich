#source('R/HallMarkAnalysis_Preamble.R')

TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV','PRAD')

for (TYPE in TYPES){

pdf(paste('../../RESULTS/SLAPenrich/ALL_radialPlots/',TYPE,'.pdf',sep=''),width = 7,height = 7.5)
#RES<-SLE.HallmarkAnalysis(TYPE=TYPE,PATH_COLLECTION = PATH_COLLECTION,HM_TABLE = RES$HallMark_Table,PATH = '../../RESULTS/SLAPenrich/PT_HM_20160718/',
#                          removeDrivers=TRUE)

load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPE,'_HM.rdata',sep=''))
tmp<- (RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
names(tmp)<-RES$Hallmark
PATT_VEC<-list(DATA=tmp,HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC,mc = 50, RADIUS = 1,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE,
                       title=TYPE)
ttmp<-tmp
un<-unique(names(tmp))
for (i in 1:length(un)){
    ttmp[which(names(tmp)==un[i])]<-rep(sum(tmp[which(names(tmp)==un[i])])/length(which(names(tmp)==un[i])),
                                        length(which(names(tmp)==un[i])))
}
PATT_VEC<-list(DATA=ttmp,HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC,maxscale=1,mc = 50, RADIUS = 1.5,TOADD = TRUE,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)

PATT_VEC<-list(DATA=rep(0,length(ttmp)),HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC,maxscale=1,mc = 50, RADIUS = 1.95,TOADD = TRUE,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)
dev.off()
}
