load('data/SLAPE.hgnc.table_20160210.Rdata')
load('data/SLAPE.20160211_MSigDB_KEGG_hugoUpdated.Rdata')
load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")

source('R/SLAPenrich.R')

#data("SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection")
#updated.hgnc.table<-SLAPE.UpdateHGNC.Table()
#save(updated.hgnc.table,file='data/SLAPE.hgnc.table_20160210.Rdata')

#PATHCOM_HUMAN<-
#    SLAPE.Check_and_fix_PathwayCollection(Pathways = PATHCOM_HUMAN,updated.hgnc.table = updated.hgnc.table)
#save(PATHCOM_HUMAN,file='data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.RData')

source('R/createLungDataset.R')

PFPw<-SLAPE.Analyse(wBEM = Dataset,PATH_COLLECTION = KEGG_PATH,
                    show_progress = TRUE,
                    NSAMPLES = 1,
                    NGENES = 1,
                    accExLength = TRUE,
                    BACKGROUNDpopulation = rownames(Dataset),
                    path_probability = 'Bernoulli',
                    GeneLenghts = GECOBLenghts)

SLAPE.write.table(PFP = PFPw,BEM = Dataset,filename = "../../RESULTS/SLAPenrich/LungDS_KEGG_enrichments.csv",
                  fdrth=5,exclcovth = 0,PATH_COLLECTION = KEGG_PATH)

SLAPE.serialPathVis(BEM = Dataset,PFP = PFPw,fdrth = 1,exCovTh = 0,PATH = "../../RESULTS/SLAPenrich/",PATH_COLLECTION = KEGG_PATH)


SLAPE.coreComponents(PFP = PFPw,BEM = Dataset,fdrth = 5,exclcovth = 0,filename = 'temp.results/core.components.pdf',PATH_COLLECTION = KEGG_PATH)
