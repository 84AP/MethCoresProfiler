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

## devtools is required
#library(devtools)
#install_github("trinker/pacman")

## obtain path
#(WD <- getwd())
#WD1=sub("/[^/]+$", "", WD)
#worDir=paste(WD1,"/testData/",sep="")
#if (!is.null(WD)) setwd(worDir)

## Set path  where you have saved MethCoresProfiler-master
#setwd("workDir") #("~/MethCoresProfiler-master/testData/")

## Set path  where you have saved MethCoresProfiler-master
#setwd("~/MethCoresProfiler-master/testData/")

args=commandArgs(trailingOnly=T)
## Load all required packages
#Packages <- c("vegan", "Hmisc", "data.table", "FunChisq", "psych", "PerformanceAnalytics", "gtools", "ggpubr", "gridExtra", "grid", "rcompanion", "plyr", "tidyr",
#              "dplyr", "entropy", "tidyverse", "plotrix", "igraph", "network", "foreach", "doParallel", "stringr", "stringi", "ggcorrplot", "corrplot",
#              "scales","ggplot2","reshape2","RelValAnalysis","RFLPtools","dendextend","ape","stringdist","reshape","rgl","factoextra","stats",
#              "mclust","ComplexHeatmap","gplots","caret","png","ggplot2","circlize","cluster","pastecs","cowplot") #"doSNOW","doMC",

#Packages %in% loadedNamespaces() # check if the packages are loaded
# [1] FALSE FALSE

#pacman::p_load(Packages, character.only = TRUE)

#Packages %in% loadedNamespaces()

## Install "rowr" Package
## Package ‘rowr’ was removed from the CRAN repository.
## Formerly available versions can be obtained from the archive by clicking on this link:https://cran.r-project.org/src/contrib/Archive/rowr/

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

source(paste(scriptsDir,"/FUN/combinations_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/summ_stat_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/combs_summary1_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/combs_summ_clust1_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/chi_square_funz.R",sep=""))

#remove all of the variables from foreach:::.foreachGlobals which is where foreach keeps all of its state
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


pvalue=pvalue1 #0.0000000001
#Type of Selection
Selection=data.frame(Selection=c(paste("Positive Selection, pvalue<=",pvalue,sep=""),"No Selection")) #paste("Negative Selection, pvalue<=",pvalue,sep="")
Selection[]=lapply(Selection, as.character)

class="Methylated"

#genes=list.files(path = ".", pattern="*.cgpos")
#maps=list.files(path = ".", pattern="*.map")

#Number of possible combinations
Complexity=data.frame(Names=c("di","tri","tetra","penta","hexa","hepta","octa","ennea","deca"), 
                      Complexity=c("2","3","4","5","6","7","8","9","10"))

#create a new folder for results
input_folder=paste(input_outputDir,"2CpGs/",sep="")
dir.create(paste(input_outputDir,"All_Combinations/",sep=""))
destination_folder=paste(input_outputDir,"All_Combinations/",sep="")
#"./",destination_folder,"/",
######################################################## non toccare ####################################################
#load clones and combination list
#Mcloni <- list.files(pattern = '.*_Methylated.clones.txt')  #carica cloni Methylated
# Load genes list
genes=list.files(path = paste(input_outputDir,"/",metaCGpos,sep=""), pattern="*.cgpos")
genes=paste(input_outputDir,"/",metaCGpos,genes,sep="")
## Load metamaps list
maps=list.files(path = paste(input_outputDir,"/",metaCGpos,sep=""), pattern="*.map")
maps=paste(input_outputDir,"/",metaCGpos,maps,sep="")
#gene="DDO1.cgpos"

nProcessor <- makeCluster(nProcessor1)
#registerdoSNOW(cl)
registerDoParallel(nProcessor)
clusterExport(nProcessor, ls())

#gene1=2
for (gene1 in 1:length(genes)) {
  #Set n processors
  # nProcessor <- makeCluster(nProcessor1)
  # #registerdoSNOW(cl)
  # registerDoParallel(nProcessor)
  #Export nCores and all files to the environment
  #clusterExport(nProcessor, ls())
  
  #Set max dim CG_pos to combination
  maxCG_pos=maxCG_pos1
  #Set pvalue
  pvalue=pvalue1 #0.0000000001
  HighComplexity=HighComplexity1
  
  
  source(paste(scriptsDir,"/FUN/combinations_funz.R",sep=""))
  source(paste(scriptsDir,"/FUN/summ_stat_funz.R",sep=""))
  source(paste(scriptsDir,"/FUN/combs_summary1_funz.R",sep=""))
  source(paste(scriptsDir,"/FUN/combs_summ_clust1_funz.R",sep=""))
  source(paste(scriptsDir,"/FUN/chi_square_funz.R",sep=""))
  
  CG_pos <- read.table(genes[gene1], quote="\"", comment.char="", stringsAsFactors=FALSE)
  gene_name=strsplit(genes[gene1], "/")[[1]]
  gene_name=strsplit(gene_name[length(gene_name)], "\\.")[[1]]
  
  #gene_name=strsplit(genes[gene1], "\\.")[[1]]
  gene=paste(gene_name[1],sep="")
  
  ## create a gene folder
  dir.create(paste(destination_folder,"/",gene,sep=""))
  gene_destination_folder=paste(destination_folder,"/",gene,sep="")
  
  map=read.delim(paste(input_outputDir,"/",metaCGpos,"/",gene,"_meta.map",sep=""), header = T,quote="\"", comment.char="", stringsAsFactors=FALSE)
  #map=read.delim(paste(input_outputDir,"/",gene,"_meta.map",sep=""), header = T,quote="\"", comment.char="", stringsAsFactors=FALSE)
  
  ## list file
  Mclust_Positive <- list.files(path= paste(input_folder,gene,"/",sep=""),pattern = '.*Positive.*._Methylated.*.txt') #carica clust_files Methylated
  Random_Mclust <- list.files(path=paste(input_folder,gene,"/",sep=""),pattern = 'Random.*.All_Links_Methylated_clust.*.txt') #carica clust_files Methylated
  Mclust=append(Mclust_Positive,Random_Mclust)
  
  #Remove Random_Positive from list
  Mclust=Mclust[!grepl("Random_Positive_", Mclust, fixed=T)]
  
  #Mclust=paste("./",input_folder,"/",Mclust,sep="")

  #tab of average methylation
  #stat=paste(input_folder,gene,"/Table_Stat_Methylation_",gene,".txt",sep="")
  Stat=read.table(paste(input_folder,gene,"/Table_Stat_Methylation_",gene,".txt",sep=""),header = T)
  colnames(Stat)=c(CG_pos[1,],"Means","sd","Samples","se")
  
  #tab of Expected freq
  Expected_Freq=as.data.frame(matrix(nrow=(dim(CG_pos)[2]-1), ncol=2))
  colnames(Expected_Freq)=c("dim_core","Freq")
  Expected_Freq$dim_core=c(2:(dim(CG_pos)[2]))
  Expected_Freq$Freq=(0.5^Expected_Freq$dim_core)
  Expected_Freq$dim_core=as.character(Expected_Freq$dim_core)
  
  #Try to upload a single file
  #h=38
  gene_Mclust=Mclust[grep(gene, Mclust, fixed=T)]
  
############# Mclust analysis 
    for(h in 1:length(gene_Mclust)) {
      
    test=data.frame(matrix(ncol=1,nrow=0))
    colnames(test)="core"
    
    #decompose the filename "i"
    a=strsplit(gene_Mclust[h], ".txt")[[1]] 
    aa=strsplit(a, "\\_")[[1]] 
    #a1=strsplit(aa, "\\_")[[1]]
    a1=aa
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
    #Selection type
    Clust_Type=paste0(a1[(length(a1)-4):(length(a1)-3)],collapse="_")
    
    #loads the "clones" file by name
    j=paste(input_folder,gene,"/",sam,"_",class,"_",gene,"_","clones.txt",sep="")
    clones=read.table(j, header=TRUE, sep="\t",stringsAsFactors=FALSE, colClasses = c(rep("character",2), rep("numeric", 5)))

    #load the file "list" of the combinations to 2
    #info = file.info(Mclust[h])
    # empty = rownames(info[info$size == 0, ])
    #if (info$size == 0) next
    list=read.table(paste(input_folder,gene,"/",gene_Mclust[h],sep=""), header=T, comment.char = "A", sep="\t")

    #if the combination list is empty, close the loop and load another sample
    if (!is.data.frame(list) || !nrow(list)) {
      final_summary=data.frame(samples=sam,core=NA,Tot=clones$tot_reads_sample,nEpialCore=-1,PercEpialCore=-1,nEpialUnici=-1,PercEpialCoreUnici=-1,PercEpialMoreFreq=-1,
                               shannon_entropy=-1,RenyiEntropy=-1,expected_Freq=-1,observed_Freq=-1,Exp_xSquared=-1,Exp_pvalue=-1,n=-1,Random_Freq=-1,
                               Random_xSquared=-1,Random_pvalue=-1,Selection=NA)
      final_summ_clust=clones[clones$n>1,]
      
      #write the dimethyl files
      write.table(x=final_summary, file=paste(gene_destination_folder,"/","ComposCore_dimethyl_",class,"_",sam,"_",Clust_Type,"_" ,gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
      write.table(x=final_summ_clust, file=paste(gene_destination_folder,"/","Tab_Epialleli_dimethyl_",class,"_", sam, "_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
      
      #write the dimethyl files
      write.table(x=final_summary, file=paste(gene_destination_folder,"/","Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
      write.table(x=final_summ_clust, file=paste(gene_destination_folder,"/","Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
   
      next 
      
       } else {
      list=list
      }
    
    # create a new column `x` with the three columns collapsed together
    cluster=data.frame(core=apply( list[ ,c("From","to") ] , 1 , paste , collapse = "-" ))
    cluster$PercEpialCore=list$observed_Freq*100
    
    #z = percentile(cl = cl)
    
    #eliminates less than average combinations
    #cluster=cluster[cluster$PercEpialCore> z,]
    
####### Dymetil
    summary2=combinations(cl=cluster,tab_Stat=Stat,gene_destination_folder)
    
    #write all cores in one file "summary_cores"
    final_summary=summary2$summary1
    
    #write all epialleli files in a single "summary_epialleli" file
    final_summ_clust=summary2$summ_clust1
    
    test=rbind(test,summary2$test2)
    
    #save files
    write.table(x=final_summary, file=paste(gene_destination_folder,"/","Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste(gene_destination_folder,"/","Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
    
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
    
    if (HighComplexity==FALSE) {
      clust=clust[1,]
    } else {
      clust=clust
    }
    
    # Check if rows ==1
    if (is.null(clust) || nrow(clust) <=2) next
    
    #z = percentile(cl = cl)
    
    #eliminates less than average combinations
    #clust=clust[clust$PercEpialCore> z,]
    
    #Generates all possible combinations
    clust1=as.data.frame(t(combn(clust$core,2, fun= NULL)))
    
    #creates an automatically joined column
    cluster=tidyr::unite_(clust1, paste(colnames(clust1), collapse="_"), colnames(clust1),sep="-")
    
    #name the combination column
    colnames(cluster) <- "core"
    
    summary2=combinations(cl=cluster,tab_Stat=Stat,gene_destination_folder)
    
    #write all cores in one file "summary_cores"
    final_summary=rbind(final_summary,summary2$summary1)
    final_summary= final_summary[!duplicated(final_summary$core),,drop=FALSE]
    
    #write all epialleli files in a single "summary_epialleli" file
    final_summ_clust=rbind(final_summ_clust,summary2$summ_clust1)
    final_summ_clust= final_summ_clust[!duplicated(final_summ_clust$core),,drop=FALSE]
    
    test=rbind(test,summary2$test2)
    
    #saves files
    write.table(x=final_summary, file=paste(gene_destination_folder,"/","Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste(gene_destination_folder,"/","Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    
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
    
    #eliminates less than average combinations
    #clust=clust[clust$PercEpialCore> z,]
    
    #Generates all possible combinations
    clust1=as.data.frame(t(combn(clust$core,2, fun= NULL)))
    
    #creates an automatically joined column
    cluster=tidyr::unite_(clust1, paste(colnames(clust1), collapse="_"), colnames(clust1),sep="-")
    
    #assign name the combination column
    colnames(cluster) <- "core"
    
    summary2=combinations(cl=cluster,tab_Stat = Stat,gene_destination_folder)
    
    #write all cores in one "summary_cores" file
    final_summary=rbind(final_summary,summary2$summary1)
    final_summary= final_summary[!duplicated(final_summary$core),,drop=FALSE]
    
    #write all epialleli files in a single "summary_epialleli" file
    final_summ_clust=rbind(final_summ_clust,summary2$summ_clust1)
    final_summ_clust= final_summ_clust[!duplicated(final_summ_clust$core),,drop=FALSE]
    
    test=rbind(test,summary2$test2)
    
    #saves the files
    write.table(x=final_summary, file=paste(gene_destination_folder,"/","Summary_ComposCore_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    write.table(x=final_summ_clust, file=paste(gene_destination_folder,"/","Summary_Tab_Epialleli_",class,"_",sam,"_",Clust_Type,"_",gene, ".txt", sep=""), row.names=F, quote = F, sep="\t")
    
    
    rm(cluster,summary2,clust,clust1)
    
    #Erase everything except
    rm(list = ls()[!ls() %in% c("RAM","nProcessor","metaCGpos","HighComplexity1","comninations","summ_stat","combs_summary1","combs_summ_clust1","chi_square","nProcessor",
                                "maxCG_pos","pvalue","Selection","class","Complexity","input_folder","destination_folder","genes","maps",
      "vertex_cutOFF1","gene_destination_folder","input_folder","destination_folder","maxCG_pos","Stat","Selection","pvalue","CG_pos","Expected_Freq","percentile","nProcessor","combs_summary1","combs_summ_clust1","genes", "maps","summ_stat","unregister","combinations",
                                "Mcloni","Mclust","gene_Mclust","Mclust","chi_square","Complexity","map","final_summary","final_summ_clust","test","class","gene1","gene","Stat",
                                "Random_Mclust","Mclust_Positive",
                                "MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                                "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                                "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])
    
    } else {
    rm(cluster,summary2,clust,clust1)
    
      #Erase everything except
      rm(list = ls()[!ls() %in% c("RAM","nProcessor","metaCGpos","HighComplexity1","comninations","summ_stat","combs_summary1","combs_summ_clust1","chi_square","nProcessor",
                                  "maxCG_pos","pvalue","Selection","class","Complexity","input_folder","destination_folder","genes","maps",
        "vertex_cutOFF1","gene_destination_folder","input_folder","destination_folder","maxCG_pos","Stat","Selection","pvalue","CG_pos","Expected_Freq","percentile","nProcessor","combs_summary1","combs_summ_clust1","genes", "maps","summ_stat","unregister","combinations",
                                 "Mcloni","Mclust","gene_Mclust","Mclust","chi_square","Complexity","map","final_summary","final_summ_clust","test","class","gene1","gene","Stat",
                                 "Random_Mclust","Mclust_Positive",
                                 "MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                                 "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                                 "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])
      
    }
   }
  
  #Erase everything except
  rm(list = ls()[!ls() %in% c("RAM","nProcessor","metaCGpos","HighComplexity1","comninations","summ_stat","combs_summary1","combs_summ_clust1","chi_square","nProcessor",
                              "maxCG_pos","pvalue","Selection","class","Complexity","input_folder","destination_folder","genes","maps",
                              "vertex_cutOFF1","MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                              "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                              "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])
  
    }
stopCluster(nProcessor)

rm(list = ls()[!ls() %in% c("RAM","metaCGpos","HighComplexity1","vertex_cutOFF1","MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                            "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                            "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])


#Error in `[.data.frame`(Mcor, , c("id_pos", SamplesID$order)) : 
#undefined columns selected


