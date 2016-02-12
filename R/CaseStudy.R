#data("SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection")
load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")
load('data/SLAPE.hgnc.table_20160210.Rdata')
load('data/SLAPE.20160211_MSigDB_KEGG_updatedHGNC.RData')
#updated.hgnc.table<-SLAPE.UpdateHGNC.Table()
#save(updated.hgnc.table,file='data/SLAPE.hgnc.table_20160210.Rdata')

#PATHCOM_HUMAN<-
#    SLAPE.Check_and_fix_PathwayCollection(Pathways = PATHCOM_HUMAN,updated.hgnc.table = updated.hgnc.table)
#save(PATHCOM_HUMAN,file='data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.RData')

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

PFPw<-SLAPE.Analyse(wBEM = Dataset,PATH_COLLECTION = KEGG_PATH,
              show_progress = TRUE,
              NSAMPLES = 2,
              NGENES = 2,
              accExLength = TRUE,
              BACKGROUNDpopulation = rownames(Dataset))

SLAPE.write.table(PFP = PFPw,BEM = Dataset,filename = "temp.results/tmp.csv",fdrth=5,exclcovth = 0,PATH_COLLECTION = KEGG_PATH)

SLAPE.serialPathVis(BEM = Dataset,PFP = PFPw,fdrth = 5,exCovTh = 0,PATH = 'temp.results/',PATH_COLLECTION = KEGG_PATH)


SLAPE.coreComponents(PFP = PFPw,BEM = Dataset,fdrth = 5,exclcovth = 0,filename = 'temp.results/core.components.pdf',PATH_COLLECTION = KEGG_PATH)

