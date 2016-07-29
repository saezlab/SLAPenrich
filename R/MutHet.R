

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

load('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/ANNOTATIONS/BOTH/R/IncludedCancerTypes.rdata')

NN<-names(nMutationsMean)
COL<-
    IncludedCancerTypeAnnotations[NN,3]
    
par(mfrow=c(1,2))
par(xpd=TRUE)
plot(nMutationsMean,nEnrichedPathways,pch=21,log='y',frame.plot = FALSE,cex.lab=1.3,
     main=paste('R =',format(cor(nMutationsMean,nEnrichedPathways),digits = 3)),ylim=c(50,200),xlim=c(0,450),
     xlab='Avg. n. of mutated genes per sample',ylab='n. of SLAPenriched pathways',cex=1.5,bg=COL,col='black',
     cex.main=1.8)
text(nMutationsMean,nEnrichedPathways,names(nMutationsMean),pos = 3,cex=1)


plot(nSamples,nEnrichedPathways,pch=21,log='y',frame.plot = FALSE,cex.lab=1.3,
     main=paste('R =',format(cor(nSamples,nEnrichedPathways),digits = 3)),ylim=c(50,200),xlim=c(0,1200),
     xlab='n. samples',ylab='',cex=1.5,bg=COL,col='black',
     cex.main=1.8)
text(nSamples,nEnrichedPathways,names(nMutationsMean),pos = 3,cex=1)


par(mfrow=c(1,1))

load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_250.rdata')
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_400.rdata')
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/BRCA_800.rdata')

PCHpatt<-c(1,1,1,16)
x<-c(250,400,800,nSamples['BRCA'])
y<-c(mean(RES_250),mean(RES_400),mean(RES_800),nEnrichedPathways['BRCA'])
sd<-c(sd(RES_250),sd(RES_400),sd(RES_800),0)
plot(x,y,type='b',lty=2,ylim=c(50,200),xlim=c(200,1200),pch=PCHpatt,cex=1.2,lwd=2,log='y',cex.lab=1.3,
     col=IncludedCancerTypeAnnotations['BRCA',3],xlab='n. samples',ylab='n. of SLAPenriched pathways')
segments(x, y-sd,x, y+sd,col=IncludedCancerTypeAnnotations['BRCA',3])
epsilon = 10
segments(x-epsilon,y-sd,x+epsilon,y-sd,col=IncludedCancerTypeAnnotations['BRCA',3])
segments(x-epsilon,y+sd,x+epsilon,y+sd,col=IncludedCancerTypeAnnotations['BRCA',3])

PCHpatt<-c(1,16)
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/COREAD_250.rdata')
lines(c(250,nSamples['COREAD']),c(mean(RES_250),nEnrichedPathways['COREAD']),type='b',lty=2,pch=PCHpatt,lwd=2,
      cex=1.2,col=IncludedCancerTypeAnnotations['COREAD',3])
x<-c(250,nSamples['COREAD'])
y<-c(mean(RES_250),nEnrichedPathways['COREAD'])
sd<-c(sd(RES_250),0)
segments(x, y-sd,x, y+sd,col=IncludedCancerTypeAnnotations['COREAD',3])
epsilon = 10
segments(x-epsilon,y-sd,x+epsilon,y-sd,col=IncludedCancerTypeAnnotations['COREAD',3])
segments(x-epsilon,y+sd,x+epsilon,y+sd,col=IncludedCancerTypeAnnotations['COREAD',3])

PCHpatt<-c(1,16)
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/GBM_250.rdata')
lines(c(250,nSamples['GBM']),c(mean(RES_250),nEnrichedPathways['GBM']),type='b',lty=2,pch=PCHpatt,lwd=2,
      cex=1.2,col=IncludedCancerTypeAnnotations['GBM',3])
x<-c(250,nSamples['GBM'])
y<-c(mean(RES_250),nEnrichedPathways['GBM'])
sd<-c(sd(RES_250),0)
segments(x, y-sd,x, y+sd,col=IncludedCancerTypeAnnotations['GBM',3])
epsilon = 10
segments(x-epsilon,y-sd,x+epsilon,y-sd,col=IncludedCancerTypeAnnotations['GBM',3])
segments(x-epsilon,y+sd,x+epsilon,y+sd,col=IncludedCancerTypeAnnotations['GBM',3])

PCHpatt<-c(1,16)
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/HNSC_250.rdata')
lines(c(250,nSamples['HNSC']),c(mean(RES_250),nEnrichedPathways['HNSC']),type='b',lty=2,pch=PCHpatt,lwd=2,
      cex=1.2,col=IncludedCancerTypeAnnotations['HNSC',3])
x<-c(250,nSamples['HNSC'])
y<-c(mean(RES_250),nEnrichedPathways['HNSC'])
sd<-c(sd(RES_250),0)
segments(x, y-sd,x, y+sd,col=IncludedCancerTypeAnnotations['HNSC',3])
epsilon = 10
segments(x-epsilon,y-sd,x+epsilon,y-sd,col=IncludedCancerTypeAnnotations['HNSC',3])
segments(x-epsilon,y+sd,x+epsilon,y+sd,col=IncludedCancerTypeAnnotations['HNSC',3])

PCHpatt<-c(1,16)
load('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/LUAD_250.rdata')
lines(c(250,nSamples['LUAD']),c(mean(RES_250),nEnrichedPathways['LUAD']),type='b',lty=2,pch=PCHpatt,lwd=2,
      cex=1.2,col=IncludedCancerTypeAnnotations['LUAD',3])
x<-c(250,nSamples['LUAD'])
y<-c(mean(RES_250),nEnrichedPathways['LUAD'])
sd<-c(sd(RES_250),0)
segments(x, y-sd,x, y+sd,col=IncludedCancerTypeAnnotations['LUAD',3])
epsilon = 10
segments(x-epsilon,y-sd,x+epsilon,y-sd,col=IncludedCancerTypeAnnotations['LUAD',3])
segments(x-epsilon,y+sd,x+epsilon,y+sd,col=IncludedCancerTypeAnnotations['LUAD',3])

legend('bottomright',legend=c('BRCA','COREAD','GBM','HNSC','LUAD'),border = FALSE,bty = 'n',
       fill =IncludedCancerTypeAnnotations[c('BRCA','COREAD','GBM','HNSC','LUAD'),3])

legend('right',legend=c('all samples','down-sampled'),bty='n',pch=c(19,1))


