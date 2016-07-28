

library(plotrix)

compute.n.mutations.per.cancertype<-function(){
    
    TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV','PRAD')
    
    n.mutations<-list()
    
    ntypes<-length(TYPES)
    
    for (i in 1:ntypes){
        TYPE<-TYPES[i]
        
        print(TYPE)
        load(paste('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/SEQUENCING/TUMOURS/R/BEMs/',TYPE,'.rdata',sep=''))
        Variants<-BEM$logic
        
        
        
        n.mutations[[i]]<-colSums(Variants)
        names(n.mutations)[i]<-TYPE
    }
    
    save(n.mutations,file='data/n.mutations.rdata')
}

load('data/n.mutations.rdata')

TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV','PRAD')

nEnrichedPathways<-vector()
nMutationsMean<-unlist(lapply(n.mutations,'mean'))
nMutationsSD<-unlist(lapply(n.mutations,'sd'))
nSamples<-unlist(lapply(n.mutations,'length'))

nEnrichedPathways<-vector()

for (i in 1:length(TYPES)){
    load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',TYPES[i],'_HM.rdata',sep=''))
    
    idx<-which(RES$pvals<0.05 & RES$FDR<5 & RES$ME>50)
    
    nEnrichedPathways[i]<-length(idx)
}

names(nEnrichedPathways)<-TYPES

par(mfrow=c(1,2))
plot(nMutationsMean,nEnrichedPathways,pch=16,log='y',frame.plot = FALSE,main=paste('R =',format(cor(nMutationsMean,nEnrichedPathways),digits = 3)),ylim=c(50,200),xlim=c(0,450))
text(nMutationsMean,nEnrichedPathways,names(nMutationsMean),pos = 3,cex=0.5)

plot(nSamples,nEnrichedPathways,pch=16,log='y',frame.plot = FALSE,main=paste('R =',format(cor(nSamples,nEnrichedPathways),digits = 3)),ylim=c(50,200),xlim=c(0,1200))
text(nSamples,nEnrichedPathways,names(nMutationsMean),pos = 3,cex=0.5)

load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_250.rdata')
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_400.rdata')
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_800.rdata')

x<-c(250,400,800,nSamples['BRCA'])
y<-c(mean(RES_250),mean(RES_400),mean(RES_800),nEnrichedPathways['BRCA'])

sd<-c(sd(RES_250),sd(RES_400),sd(RES_800),0)

plot(x,y,type='b',lty=2,ylim=c(50,200),xlim=c(0,1200),ylim=c(50,200))
segments(x, y-sd,x, y+sd)
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd)
segments(x-epsilon,y+sd,x+epsilon,y+sd)

load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/COREAD_250.rdata')
lines(c(250,nSamples['COREAD']),c(mean(RES_250),nEnrichedPathways['COREAD']),type='b',lty=2)

load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/LUAD_250.rdata')
lines(c(250,nSamples['LUAD']),c(mean(RES_250),nEnrichedPathways['LUAD']),type='b',lty=2)


