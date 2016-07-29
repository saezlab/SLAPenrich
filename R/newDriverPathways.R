source('R/HallMarkAnalysis_Preamble.R')
load('data/IntoGen_cancer_drivers_13_06_2014.rdata')
load('data/hm_colors.rdata')

type<-'SKCM'

PATH<-'../../RESULTS/SLAPenrich/newCGpathways/'

load(paste('../../../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/SEQUENCING/TUMOURS/R/BEMs/',type,'.rdata',sep=''))
Dataset<-BEM$logic
Dataset<-SLAPE.check_and_fix_gs_Dataset(Dataset,updated.hgnc.table = updated.hgnc.table)
EM<-Dataset

CD<-cancerDrivers$ActingDriver_Symbol[which(cancerDrivers$Tumor_Type==type)]

load(paste('../../RESULTS/SLAPenrich/PT_HM_20160623/',type,'_HM.rdata',sep=''))
ALLp<-RES

load(paste('../../RESULTS/SLAPenrich/PT_HM_20160718/',type,'_HM.rdata',sep=''))
NODp<-RES

id<-which(NODp$pvals<0.05 & NODp$FDR<5 & NODp$ME>50)

HMs<-sort(unique(NODp$Hallmark[id]))

for (HM in HMs){

    id<-which(NODp$Hallmark==HM & NODp$pvals<0.05 & NODp$FDR<5 & NODp$ME>50)
    
    for (i in 1:length(id)){
        IDD<-id[i]    

        genes<-PATHCOM_HUMAN$HGNC_SYMBOL[[NODp[IDD,1]]]
        
        nGenesInPath<-length(genes)
        genes<-intersect(genes,rownames(Dataset))
        
        fn<-paste(PATH,type,'_',HM,'_',IDD,'_hmap.pdf',sep='')
        
        toPlot<-matrix(c(EM[genes,]),nrow = length(genes),ncol = ncol(EM),dimnames = list(genes,colnames(EM)))
        
        toPlot<-
            SLAPE.heuristic_mut_ex_sorting(toPlot)
        
        toPlot<-toPlot[,which(colSums(toPlot)>0)]
        
        idxNC<-which(!is.element(rownames(toPlot),CD))
        
        toPlot[which(is.element(row(toPlot),idxNC) & toPlot > 0)]<-2
        
        
        NAME<-
            PATHCOM_HUMAN$PATHWAY[NODp[IDD,1]]
        
        NAME<-paste(str_trim(unlist(str_split(NAME,'//'))),collapse='\n')
        
        MAIN<-NAME
        
        MAIN<-paste(MAIN,'\n','FDR = ',format(ALLp$FDR[IDD],digits=2),'%, FDR nod = ',format(NODp$FDR[IDD],digits=2),'%',sep='')
        
        COL<-
            hm_colors[HM]
        
        COLdeb<-
            col2rgb(COL)
        
        COLdeb<-
            rgb(red=COLdeb[1],green = COLdeb[2],blue=COLdeb[3],alpha = 120,maxColorValue = 255)
        
#         pheatmap(toPlot,cluster_rows = FALSE,cluster_cols = FALSE,
#                  filename = fn,col=c('white',COLdeb,COL),
#                  legend_breaks = c(0,0.5,1,1.5,2),
#                  legend_labels=c('wt','','known HCG mut','','unknown CG mut'),
#                  main=MAIN)
        
        mutFrequencies<-100*rowSums(EM[genes,])/ncol(EM)
        HCGs<-rep(FALSE,length(mutFrequencies))
        names(HCGs)<-genes
        HCGs[intersect(names(HCGs),CD)]<-TRUE
     
        tmp<-cbind(names(HCGs),HCGs,mutFrequencies)
           
        colnames(tmp)[1]<-'node'
        
        write.table(tmp,quote=FALSE,sep='\t',row.names = FALSE,file = paste(fn,'.txt',sep=''))
    }
}

 
