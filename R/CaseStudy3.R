source('R/SLAPenrich.R')

load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")
load('data/SLAPE.hgnc.table_20160210.Rdata')
load('data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.Rdata')

source('R/createJB_dataset.R')

Dataset<-as.matrix(Dataset)

PFPw<-SLAPE.Analyse(wBEM = Dataset,PATH_COLLECTION = PATHCOM_HUMAN,
                    show_progress = TRUE,
                    NSAMPLES = 2,
                    NGENES = 2,
                    accExLength = TRUE,
                    BACKGROUNDpopulation = rownames(Dataset),
                    path_probability = 'Bernoulli',
                    GeneLenghts = GECOBLenghts)

SLAPE.write.table(PFP = PFPw,BEM = Dataset,filename = "../../RESULTS/SLAPenrich/JB_paper.csv",
                  fdrth=5,exclcovth = 50,PATH_COLLECTION = PATHCOM_HUMAN)

SLAPE.serialPathVis(BEM = Dataset,PFP = PFPw,fdrth = 5,exCovTh = 50,PATH = "../../RESULTS/SLAPenrich/JBp/",PATH_COLLECTION = PATHCOM_HUMAN)