fn<-dir('../../RESULTS/SLAPenrich/PT_HM_inigoSuggestion/')
fn<-grep('.rdata',fn,value = TRUE)

library(stringr)

tmf<-str_split(fn,'.rdata')
tmf<-unlist(tmf)[seq(1,20,2)]

fn<-paste('../../RESULTS/SLAPenrich/PT_HM_inigoSuggestion/',fn,sep='')

nf<-length(fn)


ALLnovelPaths<-rep(0,nf)
names(ALLnovelPaths)<-tmf
notLed<-rep(0,nf)
names(notLed)<-tmf

nn<-c("BRCA","GBM","COREAD","HNSC","OV","PRAD","SKCM","LUAD","THCA",'KIRC')  
nn<-sort(nn)

for (i in 1:nf){
    load(fn[i])
    ALLnovelPaths[i]<-100*length(which(c(TOTBEM)==0))/length(!is.na(c(TOTBEM)))
    notLed[i]<-100*length(which(c(TOTBEM<50)))/length(!is.na(c(TOTBEM)))

    tempMax<-matrix(NA,nrow(TOTBEM),10,dimnames = list(rownames(TOTBEM),nn))
    tempMax[,colnames(TOTBEM)]<-TOTBEM
    
    if (i == 1){
                    ttt<-tempMax
        
        }else{
            ttt<-rbind(ttt,tempMax)
            }
    }

