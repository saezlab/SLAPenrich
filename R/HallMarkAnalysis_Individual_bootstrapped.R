source('R/HallMarkAnalysis_Preamble.R')

#TYPES<-c('BRCA','COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV','PRAD')
TYPES<-c('COREAD','KIRC','LUAD','HNSC','SKCM','GBM','THCA','OV')
TYPES<-'HNSC'

load('data/n.mutations.rdata')

n.samples<-unlist(lapply(n.mutations,'length'))

for (TYPE in TYPES){
    load(paste('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/SEQUENCING/TUMOURS/R/BEMs/',TYPE,'.rdata',sep=''))
    
    Dataset<-BEM$logic
    Dataset<-SLAPE.check_and_fix_gs_Dataset(Dataset,updated.hgnc.table = updated.hgnc.table)
    
    RES_250<-SLE.HallmarkAnalysis_bootstrap(TYPE=TYPE,PATH_COLLECTION = PATH_COLLECTION,HM_TABLE = RES$HallMark_Table,PATH = '../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',
                                        NSAMPLES=250,NTRIALS = 50,Dataset = Dataset)
    
    save(RES_250,file=paste('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',TYPE,'_250.rdata',sep=''))
    
    if (TYPE=='BRCA'){
        RES_400<-SLE.HallmarkAnalysis_bootstrap(TYPE=TYPE,PATH_COLLECTION = PATH_COLLECTION,HM_TABLE = RES$HallMark_Table,PATH = '../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',
                                                NSAMPLES=400,NTRIALS = 50,Dataset = Dataset)
        save(RES_400,file=paste('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',TYPE,'_400.rdata',sep=''))
    
        RES_800<-SLE.HallmarkAnalysis_bootstrap(TYPE=TYPE,PATH_COLLECTION = PATH_COLLECTION,HM_TABLE = RES$HallMark_Table,PATH = '../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',
                                                NSAMPLES=800,NTRIALS = 50,Dataset = Dataset)
        save(RES_800,file=paste('../../RESULTS/SLAPenrich/PT_HM_20160720_bootstrapped/',TYPE,'_800.rdata',sep=''))
        
        }
}
