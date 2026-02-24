library(dplyr)
library(data.table)
library(corrplot)
library(Hmisc)
files =list.files("/data/frr/Qu_methylation/methylation_bed",pattern='*_union.bed.gz$',full=T,rec=T)
species <- c("Q.acutissima","Q.longispica","Q.rehderiana","Q.aliena","Q.aliena var. acutiserrata",
             "Q.dentata","Q.fabri","Q.griffithii","Q.robur","Q.serrata","Q.glauca","Q.glaucoides",
             "Q.mongolica","Q.chenii","Q.cocciferoides","Q.dolicholepis","Q.franchetii",
             "Q.phillyreoides","Q.variabilis")
for (i in 1:19){    #i in files
  file1 = fread(files[i])  #file1 = fread(i)
  file <-  file1 %>% select(-c(individual_num,cover_sum,methlation_ratio_mean,methlation_ratio_sum))
  list<-seq(5,ncol(file),by=7)
  meth_select <-file %>% select(unlist(list))
  species1 <- species[i]
  colnames(meth_select) <- paste(species1,1:ncol(meth_select),sep="_") 
  meth_select_re <- rcorr(as.matrix(meth_select))
  myfilename <- paste("/data/frr/Qu_methylation/analysis/plot/corrlation/", species1,".png",sep ="") 
  png(myfilename,width=4500,height=4500,res=600)
  corrplot(meth_select_re$r,  type = "upper", order = "hclust",       #移除diag = FALSE,col=COL2('RdBu', 10),col=col2(10)
           tl.pos = "d", cl.pos = "r",cl.offset=3,col=COL2('RdBu', 10),col.lim=c(0.75,1),cl.length = 6,is.corr = FALSE) 
  # add = TRUE在已有的图形上添加其他图形
  corrplot(meth_select_re$r, add = TRUE,type="lower", method = "number", order="hclust",tl.pos = "n", cl.pos = "n",
           p.mat = meth_select_re$P, sig.level = 0.05, insig = "pch",diag = FALSE,col="#2166AC" ,col.lim=c(0.75,1),is.corr = FALSE) 
  dev.off()
}

