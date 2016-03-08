variants<-read.table('../../Data/SLAPenrich/externalData/Jpaper.txt',sep='\t',row.names = 1,header=TRUE,stringsAsFactors = FALSE)


Dataset<-SLAPE.Check_and_fix_GS_Dataset(variants,updated.hgnc.table = updated.hgnc.table)
