
source('R/SLAPenrich.R')
load('../SLAPenrich/data/SLAPE.20140608_PATHCOM_HUMAN.RData')

npath<-length(PATHCOM_HUMAN$PATHWAY)

JI<-matrix(0,npath,npath)

for (i in 1:npath){
    print(i)
    g1<-
        PATHCOM_HUMAN$HGNC_SYMBOL[[i]]
    
    for (j in 1:npath){
        
        g2<-
            PATHCOM_HUMAN$HGNC_SYMBOL[[j]]
        
        JI[i,j] <-length(intersect(g1,g2))/length(union(g1,g2))
        
    }
}

JI<-JI*(upper.tri(JI)+0)

hist(c(JI[which(JI>0)]),100)

load('../SLAPenrich/data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection.Rdata')

npath<-length(PATHCOM_HUMAN$PATHWAY)

JI2<-matrix(0,npath,npath)

for (i in 1:npath){
    print(i)
    g1<-
        PATHCOM_HUMAN$HGNC_SYMBOL[[i]]
    
    for (j in 1:npath){
        
        g2<-
            PATHCOM_HUMAN$HGNC_SYMBOL[[j]]
        
        JI2[i,j] <-length(intersect(g1,g2))/length(union(g1,g2))
        
    }
}

JI2<-JI2*(upper.tri(JI2)+0)

hist(c(JI2[which(JI2>0)]),100)
