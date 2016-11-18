library(HGNChelper)
library(poibin)
library(stringr)
library(reshape)
library(png)
library(ecolitk)
source('R/SLAPenrich.R')
source('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_CODE/MISC/my.hypTest.R')

TYPE<-'SKCM'

load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPE,'_HM.rdata',sep=''))
tmp<- (RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
names(tmp)<-RES$Hallmark
PATT_VEC_ALL<-list(DATA=tmp,HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_ALL,mc = 50, RADIUS = 1,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE,
                       title=TYPE)

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
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_ALL_HM,maxscale=1,mc = 50, RADIUS = 1.5,TOADD = TRUE,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)

load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPE,'_HM.rdata',sep=''))



load(paste('../../RESULTS/SLAPenrich/PT_HM_20160718/',TYPE,'_HM.rdata',sep=''))
tmp<- (RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
names(tmp)<-RES$Hallmark
PATT_VEC_NOD<-list(DATA=tmp,HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_NOD,mc = 50, RADIUS = 1,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE,
                       title=TYPE)

ttmp<-tmp
un<-unique(names(tmp))
for (i in 1:length(un)){
    ttmp[which(names(tmp)==un[i])]<-rep(sum(tmp[which(names(tmp)==un[i])])/length(which(names(tmp)==un[i])),
                                        length(which(names(tmp)==un[i])))
}
PATT_VEC_NOD_HM<-list(DATA=ttmp,HM=rownames(RES))
SLE.radialHallmarkPlot(DATA_VEC = PATT_VEC_NOD_HM,maxscale=1,mc = 50, RADIUS = 1.5,TOADD = TRUE,SCALE = 0.5,hanna = FALSE,scaleNumb=FALSE,resort = FALSE)



pt<-unique(RES$Idx)
id<-match(pt,RES$Idx)

RES<-RES[id,]

enriched_nod<-which(RES$FDR<5 & RES$pvals<0.05 & RES$ME>50)
k<-length(enriched_nod)

x<-length(intersect(enriched_all,enriched_nod))


print(sum(PATT_VEC_ALL$DATA*PATT_VEC_NOD$DATA)/sum(PATT_VEC_ALL$DATA)*100)
print(sum(PATT_VEC_ALL$DATA*PATT_VEC_NOD$DATA)/sum(PATT_VEC_NOD$DATA)*100)

print(my.hypTest(x,k,n,N))

ALL_HM<-PATT_VEC_ALL_HM$DATA[unique(names(PATT_VEC_ALL_HM$DATA))]
NOD_HM<-PATT_VEC_NOD_HM$DATA[unique(names(PATT_VEC_NOD_HM$DATA))]
print(cor(ALL_HM,NOD_HM,method = 'spearman'))
plot(ALL_HM,NOD_HM)

