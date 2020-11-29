
#Usare se si lancia il comando da R
#setwd("~/Desktop/DRGFP_Dicembre/Tiroide_methylation/TABELLONI/RASSF1/MINUS/INV/")
#setwd("/home/Avvedimento2/DATI_LAB_VEA/Tiroide_2017/RASSF1/MINUS/out/cloni/")
setwd("/media/labveacelsius/HD4T/Chiariotti/BDNF_Chiariotti/cloni/BDNF/BDNF_Lower/")



#Usare se si lancia il comando da Terminal--> Rscript Iverti_tabella_cartella.R
#args=commandArgs(trailingOnly=T)

files=list.files(path = ".", pattern=".out")
#file="T1.out"
for (file in files) {
  sam=strsplit(x = file, split = ".out") [[1]][1]
  tab=read.table(file, quote="\"", comment.char="")
  
  #elimina CpG
  tab1=tab[,-c(1:4,18)]
  
  #inverti tabella
  ncols=dim(tab1)[2]
  tab2=rev(tab1)
   
  write.table(x = tab2, file = paste("inv", sam, ".out", sep=""), row.names=F, col.names=F, sep=" ")
}

