
library(HGNChelper)
library(poibin)
library(stringr)
library(pheatmap)

source('R/SLAPenrich.R')

load("data/SLAPE.all_genes_exonic_content_block_lengths_ensemble_20160209.RData")
load('data/SLAPE.hgnc.table_20160210.Rdata')
load('data/SLAPE.20140608_PATHCOM_HUMAN_nonredundant_intersection_hugoUpdated.Rdata')

source('R/createLungDataset.R')

load('data/caseStudy_clinicalInfos.Rdata')


SLAPE.diffSLAPE.analysis(wBEM = Dataset,contrastMatrix = CLINIC_MAT,BACKGROUNDpopulation = rownames(Dataset),SLAPE.FDRth = 5,
                         positiveCondition = 'SS_CurrentSmoker',negativeCondition = 'SS_Never',PATH_COLLECTION = PATHCOM_HUMAN)

SLAPE.diffSLAPE.analysis(wBEM = Dataset,contrastMatrix = CLINIC_MAT,BACKGROUNDpopulation = rownames(Dataset),SLAPE.FDRth = 5,
                         positiveCondition = "BAC_Type_nonMucinous",negativeCondition = "BAC_Type_Mucinous",PATH_COLLECTION = PATHCOM_HUMAN)