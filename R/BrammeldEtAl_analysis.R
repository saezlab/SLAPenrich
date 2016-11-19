source('R/SLAPenrich.R')

load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")
load('data/SLAPE.hgnc.table_20160210.Rdata')
load('data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.Rdata')
load('data/SLAPE.all_genes_exonic_lengths_ensemble_20160209.RData')
variants<-read.table('data/BrammeldEtAl_variants.txt',sep='\t',row.names = 1,header=TRUE,stringsAsFactors = FALSE)
Dataset<-
    SLAPE.check_and_fix_gs_Dataset(variants,updated.hgnc.table = updated.hgnc.table)

Dataset<-as.matrix(Dataset)

PFPw<-
    SLAPE.analyse(EM = Dataset,PATH_COLLECTION = PATHCOM_HUMAN,
                    show_progress = TRUE,
                    NSAMPLES = 2,
                    NGENES = 2,
                    accExLength = TRUE,
                    BACKGROUNDpopulation = rownames(Dataset),
                    path_probability = 'Bernoulli',
                    GeneLenghts = GECOBLenghts)

SLAPE.write.table(PFP = PFPw,EM = Dataset,filename = "../../RESULTS/SLAPenrich/JonathanPaper/JB_paper_check.csv",
                  fdrth=5,exclcovth = 50,PATH_COLLECTION = PATHCOM_HUMAN,GeneLenghts = GECOBLenghts)


SLAPE.serialPathVis(BEM = Dataset,PFP = PFPw,fdrth = 5,exCovTh = 50,PATH = "../../RESULTS/SLAPenrich/JBp/",PATH_COLLECTION = PATHCOM_HUMAN)