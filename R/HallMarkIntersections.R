par(mfrow=c(1,1))
library(pheatmap)

load('data/20160719_HallMarks_pathways_Genes.rdata')
load('data/hm_colors.rdata')

uhm<-unique(restable[,2])
upt<-unique(restable[,3])

hm_ints<-matrix(0,nrow = length(uhm),ncol = length(upt), dimnames = list(uhm,upt))

for (i in 1:nrow(restable)){
    hm_ints[restable[i,2],restable[i,3]]<-1
    }


hm_ints<-SLAPE.heuristic_mut_ex_sorting(hm_ints)*matrix(rep(1:nrow(hm_ints),ncol(hm_ints)),nrow = nrow(hm_ints),ncol=ncol(hm_ints))

pheatmap(hm_ints,show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE,legend = FALSE,col=c('white',hm_colors[rownames(hm_ints)]),border_color = 'black',main=paste(ncol(hm_ints),'Pathways'),cex=1.2)

barplot(colSums(sign(hm_ints)),names.arg = '',border = NA,col='darkgray',ylab='n. mapped Hallmarks')
barplot(rev(rowSums(sign(hm_ints))),names.arg = '',border = NA,col='darkgray',horiz = TRUE,xlab='n. mapped Pathways')


ug<-unique(restable[,5])
ug<-sort(unique(str_trim(unlist(str_split(ug,',')))))

hm_ints_genes<-matrix(0,nrow = length(uhm),ncol = length(ug), dimnames = list(uhm,ug))

for (i in 1:nrow(restable)){
    
    currentGeneSet<-restable[i,5]
    currentGeneSet<-str_trim(str_split(currentGeneSet,',')[[1]])
    
    hm_ints_genes[restable[i,2],currentGeneSet]<-1
}

hm_ints_genes<-SLAPE.heuristic_mut_ex_sorting(sign(hm_ints_genes))

hm_ints_genes<-hm_ints_genes*matrix(rep(1:nrow(hm_ints_genes),ncol(hm_ints_genes)),nrow = nrow(hm_ints_genes),ncol=ncol(hm_ints_genes))

pheatmap(hm_ints_genes,show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE,legend = FALSE,border_color = 'black',col=c('white',hm_colors[rownames(hm_ints_genes)]),main=paste(ncol(hm_ints_genes),'Genes'),cex=1.2)

barplot(colSums(sign(hm_ints_genes)),names.arg = '',border = NA,col='darkgray',ylab='n. mapped Hallmarks')
barplot(rev(rowSums(sign(hm_ints_genes))),names.arg = '',border = NA,col='darkgray',horiz = TRUE,xlab='n. mapped Pathways')


par(mfrow=c(1,2))
hist(colSums(sign(hm_ints)),col='darkgray',breaks = c(0.5,1.5,2.5,3.5,4.5),xlab='n. associated Hallmarks per Pathway',ylab='n. Pathways',main='')
hist(colSums(sign(hm_ints_genes)),col='darkgray',breaks = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab='n. associated Hallmarks per Gene',ylab='n. Genes',xlim=c(0,11),ylim=c(0,2000),main='')


restable<-restable[match(upt,restable[,3]),]
pmembmat<-matrix(0,nrow = length(upt),ncol = length(ug), dimnames = list(upt,ug))

for (i in 1:nrow(restable)){
    
    currentGeneSet<-restable[i,5]
    currentGeneSet<-str_trim(str_split(currentGeneSet,',')[[1]])
    
    pmembmat[restable[i,3],currentGeneSet]<-1
}

par(mfrow=c(1,1))
mydata_hist <- hist(1-dist(pmembmat,method = 'binary'))
barplot(mydata_hist$counts, log='y',pch=20, col="darkgray",ylab='Jaccard Index',
        main='Pathway gene set pair-wise overlaps')

axis(1)









