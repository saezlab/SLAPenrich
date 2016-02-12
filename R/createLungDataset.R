LungTSPvariants<-read.table('../../Data/SLAPenrich/externalData/supplementary_table_2.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)

GeneIds<-unique(LungTSPvariants$HUGO.Symbol)
GeneIds<-setdiff(GeneIds,'')

SampleIds<-unique(LungTSPvariants$Tumor.ID)
SampleIds<-as.character(setdiff(SampleIds,NA))

Dataset<-matrix(0,length(GeneIds),length(SampleIds),dimnames = list(GeneIds,SampleIds))


nvariants<-nrow(LungTSPvariants)
for (i in 1:nvariants){
    sampleId<-as.character(LungTSPvariants$Tumor.ID[i])
    geneId<-LungTSPvariants$HUGO.Symbol[i]
    if (geneId!='' & sampleId!='NA'){
        Dataset[geneId,sampleId]<-Dataset[geneId,sampleId]+1    
    }
}


Dataset<-SLAPE.Check_and_fix_GS_Dataset(Dataset,updated.hgnc.table = updated.hgnc.table)
