clinicInfo<-read.table('../../Data/SLAPenrich/externalData/supplementary_table_15.tsv',sep='\t',header=TRUE,stringsAsFactors = FALSE)

patientId<-
    clinicInfo$id

smokingStatus<-clinicInfo$Smoking_Status
smokingStatus[which(smokingStatus=='')]<-'SS_NotAvailable'
smokingStatus[which(smokingStatus=='F')]<-'SS_Former'
smokingStatus[which(smokingStatus=='C')]<-'SS_CurrentSmoker'
smokingStatus[which(smokingStatus=='N')]<-'SS_Never'


CLINIC_MAT<-matrix(0,nrow = length(patientId),length(unique(smokingStatus)),dimnames = list(patientId,unique(smokingStatus)))

for (i in 1:length(patientId)){
    CLINIC_MAT[as.character(patientId[i]),smokingStatus[i]]<-1
}

save(CLINIC_MAT,file='data/caseStudy_clinicalInfos.Rdata')