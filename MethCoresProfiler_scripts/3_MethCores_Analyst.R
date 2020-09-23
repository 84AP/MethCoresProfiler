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

#Keep two decimal
## Disable scientific option of R.
options("scipen"=100, "digits"=2)

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
#require("doMC")
library("progress")
library(pastecs)
require("RelValAnalysis")
library("psych")
library("cowplot")

source(paste(scriptsDir,"/FUN/compare_Samples_epialleles_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/complexity_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/strings_to_BinaryProfiles_funz.R",sep=""))
#source(paste(scriptsDir,"/FUN/epialleles_heatmap_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/cores_status_assignment_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/minimal_common_structure2_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/chi_square_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/quantile_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/MethCoresIndex_funz.R",sep=""))
#source(paste(scriptsDir,"/FUN/MethCoresIndex2_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/Assign_color_funz.R",sep=""))
#source(paste(scriptsDir,"/FUN/MethCores_tab_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/MethCores_tab2_funz.R",sep=""))
source(paste(scriptsDir,"/FUN/epialleles_heatmap2_funz.R",sep=""))


unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


#Set n processors
#detectCores()
nProcessor <- makeCluster(nProcessor1)
#registerdoSNOW(cl)
registerDoParallel(nProcessor)
clusterExport(nProcessor, ls())
#on.exit(stopCluster(nProcessor))
#stopCluster(nProcessor)

genes=list.files(path = input_outputDir, pattern="*.cgpos")
maps=list.files(path = input_outputDir, pattern="*.map")

#Set pvalue
pvalue=pvalue1 #0.0000000001

#Set max dim CG_pos to combination
maxCG_pos=maxCG_pos1

#Plot All samples togheter ("TRUE" or "FALSE")
Plot_All=Plot_All1 #"FALSE"  #or "TRUE"

#Put TRUE if you want a temporally analysis
overtime=overtime1 #"FALSE"  #or "TRUE"

saveAs=saveAs1 #pdf or png
Groups=Groups1

## if change only timing
type_Tissue=type_Tissue1 #"same" "different" 

input_folder1=paste(input_outputDir,"2CpGs/",sep="") #"input2"
input_folder2=paste(input_outputDir,"All_Combinations/",sep="") #"input3"

#create a folder
destination_folder=paste(input_outputDir,"Cores/",sep="")
dir.create(destination_folder)
### Create a Fiolder for each Group
#FOLDERS=c("CX","CB","Brain","gut","")


#Number of possible combinations
Complex=as.vector(c("di","tri","tetra","penta","hexa","hepta","octa","ennea","deca"))
HighComplexity=FALSE #TRUE or FALSE
# report weight of all epiallels contaning cores or only significant epiallels contaning cores
population_weight=population_weight1 # TRUE or FALSE
#Set cores_cutOff
cores_cutOff=c(.9,.8,.7,.6,.5,.4,.3,.2,.1)
vertex_cutOFF=0.01

class="Methylated"
class1="Methylation"

#Set which profile yoy want see
Selection="Positive"

#Options
exclude=""
select=""


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
  
  CG_pos <- read.table(paste(input_outputDir,"/",genes[gene1],sep=""), quote="\"", comment.char="", stringsAsFactors=FALSE)
  gene_name=strsplit(genes[gene1], "\\.")[[1]]
  gene=paste(gene_name[1],sep="")
  
  ## create a gene folder
  dir.create(paste(destination_folder,"/",gene,sep=""))
  gene_destination_folder=paste(destination_folder,"/",gene,sep="")
  
  map=read.delim(paste(input_outputDir,"/",gene,"_meta.map",sep=""), header = T,quote="\"", comment.char="", stringsAsFactors=FALSE)
  
  #Loading Stat
  #stat=paste(input_folder1,gene,"/","Table_Stat_Methylation_",gene,".txt",sep="")
  Stat=read.table(paste(input_folder1,gene,"/","Table_Stat_Methylation_",gene,".txt",sep=""),header = T)
  colnames(Stat)=c(CG_pos[1,],"Means","sd","Samples","se","Groups")#,"Tissue","Tissue1")
  #colnames(Stat)[dim(Stat)[2])="Tissue"
  Stat$Tissue=sub("_.*","",Stat$Groups)
  #Stat$Tissue=sub("_.*","",Stat$Groups)
  
  #Maximum complexity structures
  Complex1=Complex[c(1:(dim(CG_pos)[2]-1))]
  #Random_control=paste("Random_control_nCG",dim(CG_pos)[2],sep="")
  #Random_control=paste("Random_control_nCG",dim(CG_pos)[2],sep="")
  Random_control="Random"
  
  #loadfile CG_pos, meta.map
  #MetaMap=paste(gene,"meta.map",sep="_")
  
  #MetaMap=paste(gene1[1],"meta.map",sep="_")
  ####### Load and Check meta.map file ##########
  #load file CG_pos, meta.map
  #map=read.delim(MetaMap, stringsAsFactors=FALSE)
  map$id=paste(map$Tissue,map$Group,map$Description,map$Rep,sep="_")
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
  
  
  #Create a df as map
  tmap=map
  
  #Set Groups and remove random from list
  Groups=unique(tmap$Group)
  
  if (length(Groups)>1) {
    Groups=Groups[-(length(Groups))]
  } else {
    Groups=Groups
  }
  
  if (Plot_All=="TRUE") {
    #Groups=append(Groups,"All")
    Groups="All"
  } else {
    Groups=Groups
  }
  #Groups=Groups[-(length(Groups))]
  
  #invert columns with rows (create an array)
  map=as.data.frame(t(map))
  colnames(map)=tmap$id
  #use sample name for each column
  #colnames(map)=paste(gene1[1],map[1,],sep="_")
  
  #FOLDERS=paste(as.character(unique(Groups)),sep="")
  #FOLDERS=paste(as.character(unique(tmap$Group)),sep="") #c("CX","CB","Brain","")
  #FOLDERS=FOLDERS[-(length(FOLDERS))]
  
  colorCodes <-palette(rainbow(length(unique(tmap$Group)))) #length(Groups[-(length(Groups))])
    #colors <- c("red","white","yellow","green")
  names(colorCodes) = as.character(unique(tmap$Group))
  
  ## color for each groups
  #colorCodes <- c(CB="red", CX="green", HIPP="yellow3",
  #                CBmoe="chocolate4",CBneuro="cornflowerblue",CBastro="darkorchid3",
  #                CXmoe="chocolate4",CXneuro="cornflowerblue",CXastro="darkorchid3")
  ## select group
  
  ## use select or exclude to select / exclude a group or
 # for (iii in 1:length(FOLDERS)) {
#    if (Groups=="All") {
#      select="" #"CX"
#      exclude=""
#    } else {
#    select= FOLDERS[iii] #"CX"
#    exclude=""
#    }
    
#    if (select=="") {
      
#      destination_folder=paste(gene_destination_folder,"/",Groups,"/",sep="") 
#      dir.create(paste(gene_destination_folder,"/",Groups,"/",sep=""))
#      #destination_folder="./" 
#    } else {
      ##
#      destination_folder=paste(gene_destination_folder,"/",select,sep="")
      ##
#      dir.create(paste(gene_destination_folder,"/",select,sep=""))
#    }
    
  #Create a list of table_Epialleles files
  list <- list.files(path=paste(input_folder2,gene,"/",sep=""),pattern = '^Tab_Epialleli_.*.txt')  #carica cloni
  Random= list.files(path=paste(input_folder2,gene,"/",sep=""),pattern = '^Tab_Epialleli_.*.Random_*.All_Links.*.txt')  #carica cloni
  
  #Create a list of table_Epialleles files
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
  
  MethCores1=epialleles_heatmap2(Epi_list=Mlist,Stat,vertex_cutOFF,colorCodes,input_folder2,HighComplexity,Random,gene,gene_destination_folder,Groups)
  
  MethCores=MethCores1
  #rm(MethCores1)
  
  rm(list = ls()[!ls() %in% c("HighComplexity","input_folder1","input_folder2","type_Tissue","destination_folder","FOLDERS","colorCodes","overtime","population_weight","MethCores_tab2","saveAs","destination_folder","select","exclude","Plot_All","MethCores_tab","MethCores_tab2","Assign_color","CG_pos","Complex","Complex1","Group","tmap","Random_control","Selection","pvalue",
                              "sample_order", "maps", "genes", "list","list1","Mlist","Mplist","plist","plist1",
                              "strings_to_BinaryProfiles","compare_Samples_cores","compare_Samples_epialleles","complexity",
                              "cores_heatmap","epialleles_heatmap","epialleles_heatmap2","cores_status_assignment","minimal_common_structure2","cores_cutOff",
                              "MethCoresIndex","MethCores","MethCores_tab","maxCG_pos","quantile","Stat","class","class1","nProcessor","vertex_cutOFF",
                              "MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                              "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                              "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])  
  
 }
#}

rm(list = ls()[!ls() %in% c("MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                            "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                            "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])
