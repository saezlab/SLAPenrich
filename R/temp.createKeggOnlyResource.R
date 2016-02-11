load('data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.RData')

fc <- file('externalData/c2.cp.kegg.v5.1.symbols.txt')
mylist <- str_split(readLines(fc), "\t")
close(fc)

npath<-length(mylist)

HGNC_SYMBOL<-list()
SOURCE<-list()
PATHWAY<-list()
Ngenes<-rep(NA,npath)

for (i in 1:npath){
    pname<-mylist[[i]][1]
    
    pname<-str_split(pname,'_')[[1]]
    pname<-paste(pname[2:length(pname)],collapse=' ')
    
    
    PATHWAY[[i]]<-pname
    SOURCE[[i]]<-'KEGG'
    genes<-mylist[[i]][3:length(mylist[[i]])]
    
    HGNC_SYMBOL[[i]]<-genes
    names(HGNC_SYMBOL)[i]<-pname
    Ngenes[i]<-length(genes)
    }
names(Ngenes)<-PATHWAY

KEGG_PATH<-list(PATHWAY=PATHWAY,SOURCE=SOURCE,HGNC_SYMBOL=HGNC_SYMBOL,Ngenes=Ngenes)

KEGG_PATH<-
    SLAPE.Check_and_fix_PathwayCollection(KEGG_PATH,updated.hgnc.table = updated.hgnc.table)

save(KEGG_PATH,file = 'data/SLAPE.20160211_MSigDB_KEGG_updatedHGNC.RData')
