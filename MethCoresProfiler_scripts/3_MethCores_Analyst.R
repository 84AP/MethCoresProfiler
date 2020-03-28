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

#setwd("/media/labveacelsius/HD4T/Chiariotti/Mariella/DDO_Brain/cloni/")
#setwd("/media/labveacelsius/HD4T/Chiariotti/Mariella/Tissue_PCA/ALL/Random/")
setwd("./")

#args=commandArgs(trailingOnly=T)

#Keep two decimal
## Disabilita option scientifica di R
options("scipen"=100, "digits"=2)

#require(devtools)
require ("RFLPtools")
require ("dendextend")
require("ape")
require("stringdist")
library(reshape)
require("vegan")
library("rgl")
require("data.table")
require("factoextra")
#require("qdap")
require("stats")
require("stringi")
require("rowr")
require("entropy")
require("mclust")
require("ComplexHeatmap")
require("gplots")
require("ggpubr")
require("PerformanceAnalytics")
require("tidyr")
require("caret")
require("png")
require("ggplot2")
require("circlize")
require("tidyverse")
require("cluster")
require("stringr")
require("igraph")
require("RelValAnalysis")
library("psych")
library("cowplot")
#require("profr")
require("foreach")
require("doSNOW")
require("doParallel")
require("doMC")
library("progress")
library(pastecs)
require("RelValAnalysis")
library("psych")
library("cowplot")


#Set n processori
detectCores()
nProcessor <- makeCluster(6)
#registerdoSNOW(cl)
registerDoParallel(nProcessor)
clusterExport(nProcessor, ls())
#on.exit(stopCluster(nProcessor))
#stopCluster(nProcessor)

genes=list.files(path = ".", pattern="*.cgpos")
maps=list.files(path = ".", pattern="*.map")

#Numero di combinazioni possibili
Complex=as.vector(c("di","tri","tetra","penta","hexa","hepta","octa","ennea","deca"))

# report weight of all epiallels contaning cores or only significant epiallels contaning cores
population_weight=FALSE # TRUE or FALSE

#Set pvalue
pvalue=0.0000000001

#Set cores_cutOff
cores_cutOff=c(.9,.8,.7,.6,.5,.4,.3,.2,.1)
vertex_cutOFF=0.01

#Set max dim CG_pos to combination
maxCG_pos=8

class="Methylated"
class1="Methylation"

#Set which profile yoy want see
Selection="Positive"

#Plot All samples togheter ("TRUE" or "FALSE")
Plot_All="TRUE" #"FALSE"  #or "TRUE"

overtime="TRUE" #"FALSE"  #or "TRUE"

saveAs="pdf"#pdf or png

## select group
FOLDERS=c("CX","CB","Brain","gut","")
colorCodes <- c(CB="red", CX="green", HIPP="yellow3",
                CBmoe="chocolate4",CBneuro="cornflowerblue",CBastro="darkorchid3",
                CXmoe="chocolate4",CXneuro="cornflowerblue",CXastro="darkorchid3",
                gut="blue")
## select group
FOLDERS=c("synthetic","")
## color for each groups
colorCodes <- c(Conditionl1="red", Conditionl2="green", 
                Conditionl3="yellow3",Conditionl4="blue",
                Randomized30="cornflowerblue",Randomized10="darkorchid3")

for (iii in 1:length(FOLDERS)) {
select= FOLDERS[iii] #"CX"
exclude=""

if (select=="") {
  ##Per fare il grafico a parte
  destination_folder="./" 
} else {
  ##Per fare il grafico a parte
  destination_folder=select
  ##Crea un dataframe vuoto
  dir.create(destination_folder)
}

### Create a Fiolder for each Group
#FOLDERS=c("CX","CB","Brain","gut","")

# for (f in 1:length(FOLDERS)) {
#   
#   if (FOLDERS[f]=="") {
#     ##Per fare il grafico a parte
#     destination_folder="./" 
#   } else {
#     ##Per fare il grafico a parte
#     destination_folder=FOLDERS[f]
#     ##Crea un dataframe vuoto
#     dir.create(destination_folder)
#   }
# }
####### Sistema datafarme contenente CG_pos ##########
#gene="RandomCpGs.cgpos"
#gene1=1
for (gene1 in 1:length(genes)) {
  
  source("./Final_Script_last_26102019/compare_Samples_epialleles_funz.R")
  source("./Final_Script_last_26102019/complexity_funz.R")
  source("./Final_Script_last_26102019/strings_to_BinaryProfiles_funz.R")
  source("./Final_Script_last_26102019/epialleles_heatmap_funz.R")
  source("./Final_Script_last_26102019/cores_status_assignment_funz.R")
  source("./Final_Script_last_26102019/minimal_common_structure2_funz.R")
  source("./Final_Script_last_26102019/chi_square_funz.R")
  source("./Final_Script_last_26102019/quantile_funz.R")
  source("./Final_Script_last_26102019/MethCoresIndex_funz.R")
  source("./Final_Script_last_26102019/Assign_color_funz.R")
  source("./Final_Script_last_26102019/MethCores_tab_funz.R")
  source("./Final_Script_last_26102019/MethCores_tab2_funz.R")
  source("./Final_Script_last_26102019/epialleles_heatmap2_funz.R")
  
  
  #Carica CG_pos
  CG_pos=read.table(genes[gene1], quote="\"", comment.char="", stringsAsFactors=FALSE)
  gene_name=strsplit(genes[gene1], "\\.")[[1]]
  gene=paste(gene_name[1],sep="")
  
  #Loading Stat
  stat=paste("./Table_Stat_Methylation_",gene,".txt",sep="")
  Stat=read.table(stat,header = T)
  colnames(Stat)=c(CG_pos[1,],"Means","sd","Samples","se","Groups")#,"Tissue","Tissue1")
  #colnames(Stat)[dim(Stat)[2])="Tissue"
  Stat$Tissue=sub("_.*","",Stat$Groups)
  Stat$Tissue=sub("_.*","",Stat$Groups)
  
  #Maximum complexity structures
  Complex1=Complex[c(1:(dim(CG_pos)[2]-1))]
  #Random_control=paste("Random_control_nCG",dim(CG_pos)[2],sep="")
  #Random_control=paste("Random_control_nCG",dim(CG_pos)[2],sep="")
  Random_control="Random"
  
  ####### Sistema datafarme del meta.map ##########
  #carica file CG_pos, meta.map
  MetaMap=paste(gene,"meta.map",sep="_")
  
  #MetaMap=paste(gene1[1],"meta.map",sep="_")
  ####### Sistema datafarme del meta.map ##########
  #carica file CG_pos, meta.map
  map=read.delim(MetaMap, stringsAsFactors=FALSE)
  map$id=paste(map$Group,map$Tissue,map$Description,map$Rep,sep="_")
  map$id[dim(map)[1]]="Random"
  
  #Exclude Group Samples
  if (exclude=="") {
    map=map
    Stat=Stat
  } else {
    map=map[!grepl(exclude, map$id),]
    Stat=Stat[!grepl(exclude, Stat$Samples),]
  }
  
  #Select Group Samples
  if (select=="") {
    map=map
    Stat=Stat
  } else {
    #modify map file
    mRandom=map[grepl("Random", map$id),]
    map=map[grepl(select, map$id),]
    map=rbind(map,mRandom)
    #mofify stat file
    sRandom=Stat[grepl("Random", Stat$Samples),]
    Stat=Stat[grepl(select, Stat$Samples),]
    Stat=rbind(Stat,sRandom)
  }
  
  sample_order=as.character(map$Description)
  
  colnames(map)[1]="SamplesID"
  colnames(map)[2]="Tissue"
  colnames(map)[3]="Description"
  colnames(map)[4]="Group"
  colnames(map)[5]="Rep"
  colnames(map)[6]="id"
  
  
  #Crea a df as map
  tmap=map
  
  #Set Groups and remove random from list
  Groups=unique(tmap$Group)
  
  if (length(Groups)>1) {
    Groups=Groups[-(length(Groups))]
  } else {
    Groups=Groups
  }
  
  if (Plot_All=="TRUE") {
    Groups=append(Groups,"All")
  } else {
    Groups=Groups
  }
  #Groups=Groups[-(length(Groups))]
  
  #inverti colonne con righe (crea una matrice)
  map=as.data.frame(t(map))
  colnames(map)=tmap$id
  #usa nome campione per ogni colonna
  #colnames(map)=paste(gene1[1],map[1,],sep="_")
  
  #Crea una list dei files del gene
  list <- list.files(pattern = '^Tab_Epialleli_.*.txt')  #carica cloni
  Random= list.files(pattern = '^Tab_Epialleli_.*.Random.*.All_Links.*.txt')  #carica cloni
  
  #Crea una list dei files summary by selected
  list1=list[grep(Selection, list, fixed=T)]
  list1=list1[!grepl("Random_Positive", list1, fixed=T)]
  
  #Exclude Group Samples
  if (exclude=="") {
    list1=list1
  } else {
    list1=list1[!grepl(exclude, list1, fixed=T)]
  }
  
  #Select Group Samples
  if (select=="") {
    list1=list1
  } else {
    #lRandom=list1[grepl("Random", list1, fixed=T)]
    list1=list1[grepl(select, list1, fixed=T)]
    #list1=append(list1,lRandom)
  }
  
  # Methylated Analysis
  Mlist=list1[grep(class, list1, fixed=T)]
  #Mlist=append(Mlist,Random)
  
  ############################### HeatMap fropm Tab_Epiallels Totali
  
  MethCores1=epialleles_heatmap2(Epi_list=Mlist,Stat=Stat,vertex_cutOFF=vertex_cutOFF,colorCodes=colorCodes)
  
  MethCores=MethCores1
  #rm(MethCores1)
  
  rm(list = ls()[!ls() %in% c("FOLDERS","colorCodes","overtime","population_weight","MethCores_tab2","saveAs","destination_folder","select","exclude","Plot_All","MethCores_tab","MethCores_tab2","Assign_color","CG_pos","Complex","Complex1","Group","tmap","Random_control","Selection","pvalue",
                              "sample_order", "maps", "genes", "list","list1","Mlist","Mplist","plist","plist1",
                              "strings_to_BinaryProfiles","compare_Samples_cores","compare_Samples_epialleles","complexity",
                              "cores_heatmap","epialleles_heatmap","epialleles_heatmap2","cores_status_assignment","minimal_common_structure2","cores_cutOff",
                              "MethCoresIndex","MethCores","MethCores_tab","maxCG_pos","quantile","Stat","class","class1","nProcessor","vertex_cutOFF")])  
  
 }
}
# if (saveAs=="pdf") {
#   pdf(paste(sam, "_Tetrachoric_Correlation_", gene,".pdf", sep=""),paper = "a4")
# } else {
# }