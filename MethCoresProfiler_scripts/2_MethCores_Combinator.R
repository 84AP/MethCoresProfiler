 # MethCoresProfiler is a R-script that provides a simple method to trace and track 
 # cores shared by epiallele families in complex populations. 
 # Copyright (C) 2020 author: Antonio Pezone 
 # email: antoniopezone@gmail.com; antonio.pezone@unina.it

 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # any later version.

 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.

 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <https://www.gnu.org/licenses/>.

#setwd("/Volumes/G/cloni/")
#setwd("/media/labveacelsius/HDD2/MethCoresProfiler/AmpliMethProfilerProva/Cloni/Random_Control_Plot/")
#setwd("/media/labveacelsius/HD4T/Chiariotti/Mariella/Tissue_PCA/ALL/Random/")
setwd("./")

#args=commandArgs(trailingOnly=T)
require("RelValAnalysis")
require("dplyr")
require("plyr")
require("tidyr")
require("foreach")
require("doParallel")
require("stringr")
require("stringi")
require("entropy")
require("rcompanion")
require("data.table")
require("vegan")

source("./Final_Script_last_07102019/combinations_funz.R")
source("./Final_Script_last_07102019/summ_stat_funz.R")
source("./Final_Script_last_07102019/combs_summary1_funz.R")
source("./Final_Script_last_07102019/combs_summ_clust1_funz.R")
source("./Final_Script_last_07102019/chi_square_funz.R")

#remove all of the variables from foreach:::.foreachGlobals which is where foreach keeps all of its state
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#Set n processori
nProcessor <- makeCluster(4)
#registerdoSNOW(cl)
registerDoParallel(nProcessor)
#Esporta nCores ed all files in ambiente
clusterExport(nProcessor, ls())

#Set max dim CG_pos to combination
maxCG_pos=8
#Set pvalue
pvalue=0.0000000001
#Type of Selection
Selection=data.frame(Selection=c(paste("Positive Selection, pvalue<=",pvalue,sep=""),"No Selection")) #paste("Negative Selection, pvalue<=",pvalue,sep="")
Selection[]=lapply(Selection, as.character)

class="Methylated"

#genes=list.files(path = ".", pattern="*.cgpos")
#maps=list.files(path = ".", pattern="*.map")

#Numero di combinazioni possibili
Complexity=data.frame(Names=c("di","tri","tetra","penta","hexa","hepta","octa","ennea","deca"), 
                      Complexity=c("2","3","4","5","6","7","8","9","10"))

######################################################## non toccare ####################################################
#carica cloni e lista combinazioni
#Mcloni <- list.files(pattern = '.*_Methylated.clones.txt')  #carica cloni Methylated

Mclust_Positive <- list.files(pattern = '.*Positive.*._Methylated.*.txt') #carica clust_files Methylated

Random_Mclust <- list.files(pattern = 'Random.*.All_Links_Methylated_clust.*.txt') #carica clust_files Methylated

Mclust=append(Mclust_Positive,Random_Mclust)

genes=list.files(path = ".", pattern="*.cgpos")
maps=list.files(path = ".", pattern="*.map")
#gene="DDO1.cgpos"
#gene1=1

for (gene1 in 1:length(genes)) {

  source("./Final_Script_last_07102019/combinations_funz.R")
  source("./Final_Script_last_07102019/summ_stat_funz.R")
  source("./Final_Script_last_07102019/combs_summary1_funz.R")
  source("./Final_Script_last_07102019/combs_summ_clust1_funz.R")
  source("./Final_Script_last_07102019/chi_square_funz.R")
  
  CG_pos <- read.table(paste(genes[gene1],sep=""), quote="\"", comment.char="", stringsAsFactors=FALSE)
  gene_name=strsplit(genes[gene1], "\\.")[[1]]
  gene_names=strsplit(gene_name, "\\_")[[1]]
  gene=paste(gene_name[1],sep="")
  
  #print(genes[gene1])
#}
  #tab of map
  map=read.delim(paste(gene,"_meta.map",sep=""), stringsAsFactors=FALSE)

  #tab of average methylation
  stat=paste("./Table_Stat_Methylation_",gene,".txt",sep="")
  Stat=read.table(stat,header = T)
  colnames(Stat)=c(CG_pos[1,],"Means","sd","Samples","se")
  
  #tab of Expected freq
  Expected_Freq=as.data.frame(matrix(nrow=(dim(CG_pos)[2]-1), ncol=2))
  colnames(Expected_Freq)=c("dim_core","Freq")
  Expected_Freq$dim_core=c(2:(dim(CG_pos)[2]))
  Expected_Freq$Freq=(0.5^Expected_Freq$dim_core)
  Expected_Freq$dim_core=as.character(Expected_Freq$dim_core)
  
  #Prova per caricare un singolo file
  #h=38
  gene_Mclust=Mclust[grep(gene, Mclust, fixed=T)]
  
############# Mclust analysis 
    for(h in 1:length(gene_Mclust)) {
      
    test=data.frame(matrix(ncol=1,nrow=0))
    colnames(test)="core"
    
    #scomponi il nome del file "i"
    a=strsplit(gene_Mclust[h], "\\.")[[1]] 
    a1=strsplit(a, "\\_")[[1]]
    
   # if (length(a1)==7) {
   #   sam=paste(a1[1],"_",a1[2],sep="")
   #   #Tipo di selezione
   #   Clust_Type=paste(a1[3],"_",a1[4],sep="")
   # } else if (length(a1)==8) {
   #    sam=paste(a1[1],"_",a1[2],"_",a1[3],sep="")
   #    #Tipo di selezione
   #    Clust_Type=paste(a1[4],"_",a1[5],sep="")
   # } 
    
    sam=paste0(a1[1:(length(a1)-5)],collapse="_")
    #Tipo di selezione
    Clust_Type=paste0(a1[(length(a1)-4):(length(a1)-3)],collapse="_")
    
    #carica il file "clones" in base al nome
    j=paste(sam,class,gene,"clones.txt",sep="_")
    clones=read.table(j, header=TRUE, sep="\t",stringsAsFactors=FALSE, colClasses = c(rep("character",2), rep("numeric", 5)))

    #carica il file "lista" delle combinazioni a 2
    #info = file.info(Mclust[h])
    # empty = rownames(info[info$size == 0, ])
    #if (info$size == 0) next
    #carica il file "lista" delle combinazioni a 2
    list=read.table(gene_Mclust[h], header=T, comment.char = "A", sep="\t")

    #se la lista delle combinazioni è vuota chiudi il loop e carica altro campione
    if (!is.data.frame(list) || !nrow(list)) {
      final_summary=data.frame(samples=sam,core=NA,Tot=clones$tot_reads_sample,nEpialCore=-1,PercEpialCore=-1,nEpialUnici=-1,PercEpialCoreUnici=-1,PercEpialMoreFreq=-1,
                               shannon_entropy=-1,RenyiEntropy=-1,expected_Freq=-1,observed_Freq=-1,Exp_xSquared=-1,Exp_pvalue=-1,n=-1,Random_Freq=-1,
                               Random_xSquared=-1,Random_pvalue=-1,Selection=NA)
      final_summ_clust=clones[clones$n>1,]
      
      #scrivi i files del dimetile
      write.table(x=final_summary, file=paste("./ComposCore_dimethyl_",class,"_",sam,"_",Clust_Type,"_" ,gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
      write.table(x=final_summ_clust, file=paste("./Tab_Epialleli_dimethyl_",class,"_", sam, "_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
      
      #scrivi i files del dimetile
      write.table(x=final_summary, file=paste("./Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
      write.table(x=final_summ_clust, file=paste("./Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
      
      next
    } else {
      list=list
      }
    
    # create a new column `x` with the three columns collapsed together
    cluster=data.frame(core=apply( list[ ,c("From","to") ] , 1 , paste , collapse = "-" ))
    cluster$PercEpialCore=list$observed_Freq*100
    
    #z = percentile(cl = cl)
    
    #elimina combinazioni minori della media
    #cluster=cluster[cluster$PercEpialCore> z,]
    
####### Dymetil
    summary2=combinations(cl=cluster,tab_Stat=Stat)
    
    #scrivi tutti i cores in un unico file "summary_cores"
    final_summary=summary2$summary1
    
    #scrivi tutti i files epialleli in un unico file "summary_epialleli"
    final_summ_clust=summary2$summ_clust1
    
    test=rbind(test,summary2$test2)
    
    #salva i files
    write.table(x=final_summary, file=paste("./Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste("./Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
    
    summary1=summary2$summary1

    rm(cluster,summary2)
    
######## TetraMethyl   
    #clust=summary1[summary1$Selection %in% Selection$Selection,]
    if (stri_detect_fixed(sam,"Random")) {
    clust=summary1[summary1$Selection %in% Selection$Selection[2],]
    } else {
    clust=summary1[summary1$Selection %in% Selection$Selection[1],]
    }
    
    rm(summary1)
    
    # Check if rows ==1
    if (is.null(clust) || nrow(clust) <=2) next
    
    #z = percentile(cl = cl)
    
    #elimina combinazioni minori della media
    #clust=clust[clust$PercEpialCore> z,]
    
    #Genera tutte le combinazioni possibili
    clust1=as.data.frame(t(combn(clust$core,2, fun= NULL)))
    
    #crea una colonna unita automaticamente
    cluster=tidyr::unite_(clust1, paste(colnames(clust1), collapse="_"), colnames(clust1),sep="-")
    
    #assegna nome alla colonna combinazione
    colnames(cluster) <- "core"
    
    summary2=combinations(cl=cluster,tab_Stat=Stat)
    
    #scrivi tutti i cores in un unico file "summary_cores"
    final_summary=rbind(final_summary,summary2$summary1)
    final_summary= final_summary[!duplicated(final_summary$core),,drop=FALSE]
    
    #scrivi tutti i files epialleli in un unico file "summary_epialleli"
    final_summ_clust=rbind(final_summ_clust,summary2$summ_clust1)
    final_summ_clust= final_summ_clust[!duplicated(final_summ_clust$core),,drop=FALSE]
    
    test=rbind(test,summary2$test2)
    
    #salva i files
    write.table(x=final_summary, file=paste("./Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste("./Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    
    summary1=summary2$summary1
    
    rm(cluster,summary2,clust,clust1)
    
    # Check dim CG_Pos
    if (dim(CG_pos)[2]<maxCG_pos) {

    ######## OctoMethyl   
    #clust=summary1[summary1$Selection %in% Selection$Selection,]
      if (stri_detect_fixed(sam,"Random")) {
        clust=summary1[summary1$Selection %in% Selection$Selection[2],]
      } else {
        clust=summary1[summary1$Selection %in% Selection$Selection[1],]
      }
    
    rm(summary1)
    
    # Check if rows ==1
    if (is.null(clust) || nrow(clust) <=2) next
    
    #z = percentile(cl = cl)
    
    #elimina combinazioni minori della media
    #clust=clust[clust$PercEpialCore> z,]
    
    #Genera tutte le combinazioni possibili
    clust1=as.data.frame(t(combn(clust$core,2, fun= NULL)))
    
    #crea una colonna unita automaticamente
    cluster=tidyr::unite_(clust1, paste(colnames(clust1), collapse="_"), colnames(clust1),sep="-")
    
    #assegna nome alla colonna combinazione
    colnames(cluster) <- "core"
    
    summary2=combinations(cl=cluster,tab_Stat = Stat)
    
    #scrivi tutti i cores in un unico file "summary_cores"
    final_summary=rbind(final_summary,summary2$summary1)
    final_summary= final_summary[!duplicated(final_summary$core),,drop=FALSE]
    
    #scrivi tutti i files epialleli in un unico file "summary_epialleli"
    final_summ_clust=rbind(final_summ_clust,summary2$summ_clust1)
    final_summ_clust= final_summ_clust[!duplicated(final_summ_clust$core),,drop=FALSE]
    
    test=rbind(test,summary2$test2)
    
    #salva i files
    write.table(x=final_summary, file=paste("./Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste("./Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    
    
    rm(cluster,summary2,clust,clust1)
    
    #Cancella tutto tranne
    rm(list = ls()[!ls() %in% c("maxCG_pos","Stat","Selection","pvalue","CG_pos","Expected_Freq","percentile","nProcessor","combs_summary1","combs_summ_clust1","genes", "maps","summ_stat","unregister","combinations",
                                "Mcloni","Mclust","gene_Mclust","Mclust","chi_square","Complexity","map","final_summary","final_summ_clust","test","class","gene1","gene","Stat",
                                "Random_Mclust","Mclust_Positive")])
    
    } else {
    rm(cluster,summary2,clust,clust1)
    
      #Cancella tutto tranne
      rm(list = ls()[!ls() %in% c("maxCG_pos","Stat","Selection","pvalue","CG_pos","Expected_Freq","percentile","nProcessor","combs_summary1","combs_summ_clust1","genes", "maps","summ_stat","unregister","combinations",
                                 "Mcloni","Mclust","gene_Mclust","Mclust","chi_square","Complexity","map","final_summary","final_summ_clust","test","class","gene1","gene","Stat",
                                 "Random_Mclust","Mclust_Positive")])
      
    }
   }
  
  #Cancella tutto tranne
  rm(list = ls()[!ls() %in% c("maxCG_pos","Selection","pvalue","percentile","nProcessor","combs_summary1","combs_summ_clust1","genes", "maps","summ_stat","unregister","combinations",
                              "Mcloni","Mclust","Mclust","chi_square","Complexity","map","final_summary","final_summ_clust","class",#"gene1","gene_Mclust",
                              "Random_Mclust","Mclust_Positive")])
  
}



