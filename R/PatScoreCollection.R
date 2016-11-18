library(stringr)

filename<-'../../Data/SLAPenrich/externalData/c2.cp.v5.2.symbols.txt'

con=file(filename,open="r")
line=readLines(con) 
long=length(line)

PATHWAY<-vector()
URL<-vector()
HGNC_SYMBOL<-list()
Ngenes<-vector()
for (i in 1:long){
    print(i)
    linn=line[i]
    
    tokens<-str_split(linn,'\t')
    PATHWAY[i]<-tokens[[1]][1]
    URL[i]<-tokens[[1]][2]
    sl<-length(tokens[[1]])
    HGNC_SYMBOL[[i]]<-tokens[[1]][3:sl]
    Ngenes[i]<-length(HGNC_SYMBOL[[i]])
}

BACKGROUND<-unique(unlist(HGNC_SYMBOL))
close(con)

PATHCOM_HUMAN<-list()

PATHCOM_HUMAN$PATHWAY<-
    PATHWAY

PATHCOM_HUMAN$HGNC_SYMBOL<-
    HGNC_SYMBOL

PATHCOM_HUMAN$Ngenes<-
    Ngenes

PATHCOM_HUMAN$backGround<-
    BACKGROUND

PATHCOM_HUMAN$URL<-
    URL


save(PATHCOM_HUMAN,file='data/PATHSCORE_MsigDB.Rdata')

