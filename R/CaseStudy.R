LungTSPvariants<-read.table('externalData/supplementary_table_2.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)

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

Dataset<-SLAPE.Check_and_fix_GS_Dataset(Dataset)

#library(SLAPenrich)
data("SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection")
load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")

PATHCOM_HUMAN<-
    SLAPE.Check_and_fix_PathwayCollection(Pathways = PATHCOM_HUMAN)

PFPw<-SLAPE.Analyse(wBEM = Dataset,
              show_progress = TRUE,
              NSAMPLES = 2,
              NGENES = 2,
              accExLength = TRUE,
              BACKGROUNDpopulation = rownames(Dataset))

SLAPE.write.table(PFPw,Dataset,"tmp.txt",fdrth=5,0)
