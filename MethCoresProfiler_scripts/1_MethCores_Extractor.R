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
library(devtools)
install_github("trinker/pacman")

## obtain path
(WD <- getwd())
WD1=sub("/[^/]+$", "", WD)
worDir=paste(WD1,"/testData/",sep="")
if (!is.null(WD)) setwd(worDir)

## Set path  where you have saved MethCoresProfiler-master
#setwd("~/MethCoresProfiler-master/testData/")

## Load all required packages
Packages <- c("vegan", "Hmisc", "data.table", "FunChisq", "psych", "PerformanceAnalytics", "gtools", "ggpubr", "gridExtra", "grid", "rcompanion", "plyr", "tidyr",
              "dplyr", "entropy", "tidyverse", "plotrix", "igraph", "network", "foreach", "doParallel", "Hmisc", "stringr", "stringi", "ggcorrplot", "corrplot",
              "scales","ggplot2","reshape2","RelValAnalysis")
Packages %in% loadedNamespaces() # check if the packages are loaded
# [1] FALSE FALSE

pacman::p_load(Packages, character.only = TRUE)

Packages %in% loadedNamespaces()

#args=commandArgs(trailingOnly=T)

require("vegan")
require("Hmisc")
require("data.table")
require("FunChisq")
require("psych")
require("PerformanceAnalytics")
require("gtools")
require("ggpubr")
require("gridExtra")
require("grid")
require("rcompanion")
require("plyr")
require("tidyr")
require("dplyr")
require("entropy")
require("tidyverse")
require("plotrix")
require("igraph")
require("network")
require("foreach")
require("doParallel")
require("stringr")
require("stringi")
require("ggcorrplot")
require("corrplot")
require("scales")
require("ggplot2")
require("reshape2")
require("RelValAnalysis")

## Disable scientific option of R.
options("scipen"=100, "digits"=18)

#Set n processors
nProcessor <- makeCluster(10)
#registerdoSNOW(cl)
registerDoParallel(nProcessor)
on.exit(stopCluster(nProcessor))
#stopCluster()

genes=list.files(path = ".", pattern="*.cgpos")
maps=list.files(path = ".", pattern="*.map")

############### Set parameters
#Bootstrapping number to do
b4="default" #"default" or "fix a number, example:10000"
y=100
#Set pvalue
pvalue=0.0000000001
#Set average reads or smallest n reads
Reads="average" ## "smallest"
Remove="no"  #yes
min=300
#Set max dim CG_pos to combination
maxCG_pos=8
class="Methylated"
#Plot All samples togheter ("TRUE" or "FALSE")
Plot_All="TRUE"#FALSE"  #or "TRUE"
#make a correlation between 1CpG and 2CpG == Meth and unMeth
Meth_unMeth_corr="NO"#"YES" #or "NO
saveAs="png"#pdf or png
Modify_mono="FALSE" #"TRUE" or "FALSE" #if TRUE modify monomethyl
splitForGene="FALSE" #"splitForGene",

#create a new folder for results
dir.create("2CpGs")
destination_folder="2CpGs" #"destination_folder",
###################################### NON TOCCARE #############################################################################
for (gene1 in 1:length(genes)) {
  #gene1="p14R2.cgpos"
  #gene1=1
  
  source("../MethCoresProfiler_scripts/FUN/tetrachoric_corr_funz.R")
  source("../MethCoresProfiler_scripts/FUN/co_Meth_occurence_funz.R")
  source("../MethCoresProfiler_scripts/FUN/Mtab_clones_funz.R")
  source("../MethCoresProfiler_scripts/FUN/tab_meth_funz.R")
  source("../MethCoresProfiler_scripts/FUN/cores_extraction_funz.R")
  source("../MethCoresProfiler_scripts/FUN/make_graph_funz.R")
  source("../MethCoresProfiler_scripts/FUN/summ_stat_funz.R")
  source("../MethCoresProfiler_scripts/FUN/radian.rescale_funz.R")
  source("../MethCoresProfiler_scripts/FUN/flattenCorrMatrix_funz.R")
  source("../MethCoresProfiler_scripts/FUN/quantile_funz.R")
  source("../MethCoresProfiler_scripts/FUN/chi_square_funz.R")
  source("../MethCoresProfiler_scripts/FUN/co_Meth_unMeth_occurence_funz.R")
  source("../MethCoresProfiler_scripts/FUN/co_unMeth_Meth_occurence_funz.R")
  source("../MethCoresProfiler_scripts/FUN/correlation_plot_funz.R")
  source("../MethCoresProfiler_scripts/FUN/dectobin_funz.R")
  
  CG_pos <- read.table(paste(genes[gene1],sep=""), quote="\"", comment.char="", stringsAsFactors=FALSE)
  gene_name=strsplit(genes[gene1], "\\.")[[1]]
  gene=paste(gene_name[1],sep="")
  
  map=read.delim(paste(gene,"_meta.map",sep=""), header = T,quote="\"", comment.char="", stringsAsFactors=FALSE)
  
  # create a df with Expected Freq according CpGs
  Expected_Freq=as.data.frame(matrix(nrow=(dim(CG_pos)[2]-1), ncol=2))
  colnames(Expected_Freq)=c("dim_core","Freq")
  Expected_Freq$dim_core=c(2:(dim(CG_pos)[2]))
  Expected_Freq$Freq=(0.5^Expected_Freq$dim_core)
  Expected_Freq$dim_core=as.character(Expected_Freq$dim_core)
  
  ############### Meth_Cores_Extractor ##################################
    #Create a random control with all possibile epialleles
    #n <- as.numeric(dim(CG_pos)[2])
    #l <- rep(list(0:1), n)
    #Random=expand.grid(l)
    #colnames(Random)=CG_pos[1,]
    
  #Create a list of files.out
  outs=list.files(path = "./out", pattern=".out")
  
  if (splitForGene=="TRUE") {
    outs1=outs[grep(gene, outs, fixed=T)]
  } else {
    outs1=outs
  }
  
  outs1=paste("./out/",outs1,sep="")
  
  #Load a single file to do a test
  tab=data.frame(Tot_Reads=-1,SampleID=NA)
    
    ####################### rarefaction 
    for (out in 1:length(outs1)) { 
      #Load a file.out
      experiment <- read.table(outs1[out], quote="\"", stringsAsFactors = F)
      
      if (Modify_mono=="YES") {
        nomi_cg=as.character(CG_pos[1,])
        experiment$nMeth=rowSums( experiment[,nomi_cg]==1,na.rm=T)
        #Replace monomethyl with unmethylated
        for (i in 1:nrow( experiment)) {
          if ( experiment$nMeth[i]==1) {
            experiment[i,1:dim(CG_pos)[2]][ experiment[i,1:dim(CG_pos)[2]] == 1] <-0
          } else {
            experiment[i,1:dim(CG_pos)[2]]= experiment[i,1:dim(CG_pos)[2]]
          } 
        }
        experiment$nMeth=NULL
      } else {
        experiment= experiment
      }
      
      #Calculate the number of rows (total profiles) of the loaded sample
      #min_righe=data.frame(min(c(dim(experiment)[1])))
      tab1=data.frame(min(c(dim(experiment)[1])))
      #  min_righe=dim(experiment)[1])
      colnames(tab1)="Tot_Reads"
      #Split filename
      
      a=strsplit(outs1[out], "/")[[1]]
      a1=strsplit(a[3], "\\.")[[1]]
      #save filename as as sam
      
      sam=a1[1]
      
      tab1$SampleID=sam
      tab=rbind(tab,tab1)
    }
    
    tab=tab[!is.na(tab$SampleID),]
    
    #delete lines containing alphanumeric characters
    b=subset(tab, grepl('^\\d+$', tab$Tot_Reads))
    
    #Convert characters to numbers
    b=data.frame(b[,'Tot_Reads'] <- as.numeric(as.character(b[,'Tot_Reads'])))
    colnames(b)="Tot_Reads"
    
    #Order value
    b=data.frame(b[order(b),])
    
    # calculate average of reads or take a smallest value by parameters
    if (Reads=="average") {
      #calculates the average of the dataframe "b"
      bb4=round(as.numeric(colMeans(b, na.rm = FALSE, dims = 1)),digits = 0)
    } else { 
      bb4=round(as.numeric(b[1,1]),digits = 0)
    }
    
    if (b4=="default") {
      b4=bb4
    } else {
      b4=b4
    }
    
    # Remove smallest samples by parameters
    if (Remove=="yes") {
      #calculates the average of the dataframe "b"
      removeSample=data.frame(tab[tab$Tot_Reads <min,])$SampleID
      outs1=outs1[!grepl(removeSample, outs1)]
    } else { 
      outs1= outs1
    }
    
    #Add Random to tab
    n=dim(CG_pos)[2]
    Random_control=data.frame(Tot_Reads=(2^n),SampleID=(paste("Random_control_nCG",(dim(CG_pos)[2]),sep="")))
    rownames(Random_control)="Random"
    tab=rbind(tab,Random_control)
    tab$SamplingSize=y
    tab$Analyzed_Reads=b4
    
    tab$Tissue=map$Tissue[match(tab$SampleID, map$SampleID)]
    tab$Description=map$Description[match(tab$SampleID, map$SampleID)]
    tab$Group=map$Group[match(tab$SampleID, map$SampleID)]
    tab$Gene=gene
    tab$nCpGs=as.numeric(dim(CG_pos)[2])
    
    tab["Random","Gene"]=paste("All combinations of ",(dim(CG_pos)[2])," CpGs",sep="")
    tab=tab[ order(match(tab$SampleID, map$SampleID)), ]
    
    write.table(x = tab, file = paste("./",destination_folder,"/","Table_Reads_",gene,".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    #Delete everything and leave only the.out files
    rm(list = ls()[!ls() %in% c("splitForGene","destination_folder","Modify_mono","Meth_unMeth_corr","dectobin","chi_square","flattenCorrMatrix","radian.rescale","quantile","nProcessor","minimal_common_structure","cores_status_assignment",
                                "cores_heatmap","epialleles_heatmap","strings_to_BinaryProfiles","complexity","compare_Samples_epialleles","combs_summ_clust1",
                                "combs_summary1","combinations","summ_stat","make_graph","cores_extraction","tab_meth","Mtab_clones","co_Meth_occurence","tetrachoric_corr",
                                "nProcessor","genes","maps","Complexity","y","pvalue","Reads","Remove","min","maxCG_pos","Selection","class","Selection1",
                                "co_Meth_unMeth_occurence","co_unMeth_Meth_occurence","correlation_plot",
                                "CG_pos","gene","map","Expected_Freq","Plot_All",
                                "b4","Random","outs1","saveAs")])
    
    ########################## tab_clone_function 
    #Create a df for
    MShannon_entropy_tab=data.frame(matrix(ncol = (dim(map)[1]), nrow = 0))
    colnames(MShannon_entropy_tab)=map$X.SampleID
    
    MRenyi_entropy_tab=data.frame(matrix(ncol = (dim(map)[1]), nrow = 0))
    colnames(MRenyi_entropy_tab)=map$X.SampleID
    
    Mcor=data.frame(id_pos=NA)
    
    Mtaxa_tab=data.frame(Samples=NA,n=-1,Freq=-1,sd=-1) #n_readsGroupForSample=-1,
    
    Mepiallelic_tab=data.frame(clone=NA,id_pos=NA)
    
    ner2=as.data.frame(matrix(ncol = dim(CG_pos)[2]))
    colnames(ner2)=CG_pos[1,] 
    ner2$Means=-1
    ner2$sd=-1
    ner2$Samples=NA
    
    outs1=append(outs1,"Random")
    
    for (out in 1:length(outs1)) {
      if (outs1[out]=="Random") {
        #sam=paste("Random_control_nCG",as.numeric(dim(CG_pos)[2]),sep="")
        sam="Random"
        #n=dim(CG_pos)[2]
       # Random=sample(2^n,b4,replace=FALSE) - 1
       # Random <- Random + 2^(n)
       # names(Random) <- paste0(1:length(Random))
        #experiment <- Random %>% 
                      #  purrr::map(function(x){ dectobin(x)}) %>% 
                       # do.call(rbind,.) %>% as.data.frame() %>% select(-1)
        #colnames(experiment)=CG_pos[1,]
      } else {
        sam1=strsplit(outs1[out], "\\/")[[1]]
        sam2=strsplit(x = sam1[3], split = ".out") [[1]][1]
        sample=map[map$SampleID %in% sam2,]
        sam=paste(sample$Group,sample$Tissue,sample$Description,sample$Rep,sep="_")
        #experiment <- read.table(outs1[out], quote="\"", stringsAsFactors = F)
      }
      
      #Load tab di grandezza b4
      #tab=experiment[sample(x = dim(experiment)[1], size = b4, replace = T),]
      #tab[] <- lapply(tab[], as.numeric)
      #rownames(tab)=1:dim(tab)[1]
      #colnames(tab)=CG_pos[1,]
      
      # For Methylated CpGs
      #Create an empty lists
      tab1 <- list()
      tcs <- list()
      tetrachoric <- list()
      mclust <- list()
      unMeth_Methclust <- list()
      Meth_unMethclust <- list()
      
      #Average Methylation tab
      tab1=foreach(i=1:y,.packages=c("stringr","psych","dplyr")) %dopar% {
        tab1[[i]] = tab_meth(sam = sam, CG_pos = CG_pos,b4=b4)
      }
      #Epiaplotypes (Methylated profiles) tab
      tcs=foreach(i=1:y,.packages=c("stringr","psych","dplyr")) %dopar% {
        tcs[[i]] = Mtab_clones(sam = sam,CG_pos = CG_pos, b4=b4)
      }
      #tetrachoric tab
      tetrachoric=foreach(i=1:y,.packages=c("stringr","psych","dplyr")) %dopar% {
        tetrachoric[[i]]= tetrachoric_corr(sam = sam, CG_pos = CG_pos,b4=b4)
      }
      #co_Meth tab
      mclust=foreach(i=1:y,.packages=c("RColorBrewer","corrplot","dplyr")) %dopar% {#.options.snow=opts
        mclust[[i]]= co_Meth_occurence(sam = sam,CG_pos = CG_pos,b4=b4)
      }
      
      if (Meth_unMeth_corr=="YES") {
      #co_Meth_unMeth tab
      Meth_unMethclust=foreach(i=1:y,.packages=c("RColorBrewer","corrplot","dplyr")) %dopar% {#.options.snow=opts
        Meth_unMethclust[[i]]= co_Meth_unMeth_occurence(sam = sam,CG_pos = CG_pos,b4=b4)
      }
      #co_unMeth_Meth tab
      unMeth_Methclust=foreach(i=1:y,.packages=c("RColorBrewer","corrplot","dplyr")) %dopar% {#.options.snow=opts
        unMeth_Methclust[[i]]= co_unMeth_Meth_occurence(sam = sam,CG_pos = CG_pos,b4=b4)
      }
      
      ################# Meth_unMeth_correlation
      MethCorr=correlation_plot(tclust=mclust)
      Meth_unMethCorr=correlation_plot(tclust=Meth_unMethclust)
      unMeth_MethCorr=correlation_plot(tclust=unMeth_Methclust)
      } else {
        ################# Meth_unMeth_correlation
        MethCorr=correlation_plot(tclust=mclust)
      }
      
      
      ############## Calcolate average for tab_clones ##############
      ner=as.data.frame(matrix(ncol = dim(CG_pos)[2]))
      colnames(ner)=CG_pos[1,] 
      #i=1
      
      for (i in 1:length(tab1)){
        ner= rbind(ner,tab1[[i]])
      }
      ner = ner[complete.cases(ner), ]
      ner1=data.frame(t(colMeans(ner)))
      ner1sd=data.frame(t(apply(ner, 2, sd)))
      colnames(ner1)=CG_pos[1,]
      colnames(ner1sd)=CG_pos[1,]
      ner1f=rbind(ner1,ner1sd)
      ner1$Means=rowMeans(ner1)
      ner1$sd=sd(ner1)
      ner1$Samples=sam
      ner2=rbind(ner1,ner2)
      
      ######### CpGs Correlation ############
      #colnames(tab)=CG_pos[1,]
      #png(paste(sam,"_CpGs_Correlation_", gene,".png", sep=""), width=5*300,height=5*300, res = 300,pointsize=8)
      #chart.Correlation(tab[1:length(tab)],histogram = TRUE,method ="pearson")
      #mtext(paste(sam,"_CpGs_Correlation_", gene,sep=""), outer=TRUE,  cex=1, line=-1)
      #dev.off()
      
      ######### Tetrachoric CpGs Correlation ############
      tetra= Reduce("+",tetrachoric)/length(tetrachoric)

      if (saveAs=="pdf") {
      pdf(paste("./",destination_folder,"/",sam, "_Tetrachoric_Correlation_", gene,".pdf", sep=""),paper = "a4")
      } else {
      png(paste("./",destination_folder,"/",sam, "_Tetrachoric_Correlation_", gene,".png", sep=""), width=5*300,height=5*300, res = 300,pointsize=8)
      }
      corrplot(tetra, type="upper", order="original", sig.level = 0.01, insig = "blank", tl.srt = 100, title=paste(sam, "_Tetrachoric_Correlation_", gene,sep=""),  mar=c(0,0,1,0))
      dev.off()
      
      ###### Methylated Cores extraction
      meth_list=cores_extraction(ttcs=tcs,tclust=mclust)
      
      mer=meth_list$mer
      mer2=meth_list$mer2
      g=make_graph(tmer=mer,tmer2=mer2,sam=sam,summ_stat=summ_stat)
      
      MShannon_entropy_tab=rbind(MShannon_entropy_tab,meth_list$Shannon_entropy_tab)
      MRenyi_entropy_tab=rbind(MRenyi_entropy_tab,meth_list$Renyi_entropy_tab)
      Mcor=merge(Mcor,meth_list$cor,by="id_pos",all = T)
      Mtaxa_tab=rbind(Mtaxa_tab,meth_list$taxa_tab)
      Mepiallelic_tab=merge(Mepiallelic_tab,meth_list$epiallelic_tab,by=c("clone","id_pos"),all = T)
      
      
      #delete meth_list
      rm(meth_list)
      
      #Delete everything and leave only the.out files
      rm(list = ls()[!ls() %in% c("splitForGene","destination_folder","Modify_mono","Meth_unMeth_corr","dectobin","chi_square","flattenCorrMatrix","radian.rescale","quantile","nProcessor","minimal_common_structure","cores_status_assignment",
                                  "cores_heatmap","epialleles_heatmap","strings_to_BinaryProfiles","complexity","compare_Samples_epialleles","combs_summ_clust1",
                                  "combs_summary1","combinations","summ_stat","make_graph","cores_extraction","tab_meth","Mtab_clones","co_Meth_occurence","tetrachoric_corr",
                                  "nProcessor","genes","maps","Complexity","y","pvalue","Reads","Remove","min","maxCG_pos","Selection","class","Selection1", 
                                  "CG_pos","gene","map","Expected_Freq","Plot_All","saveAs",
                                  "co_Meth_unMeth_occurence","co_unMeth_Meth_occurence","correlation_plot",
                                  "b4","Random","ner2","Mcor","MShannon_entropy_tab","Mtaxa_tab","MRenyi_entropy_tab","Mepiallelic_tab","outs1")])
    }
    
    #################################### BarPlot sample mean methylation with SE 
    #Per fare il grafico a parte
    #ner2=read.table("./Table_Stat_Methylation_DDO1.txt",header = T, check.names = F)
    #Mtaxa_tab=read.table("./Table_Stat_Taxa_DDO1.txt",header = T, check.names = F)
    #Mcor=read.table("./Table_Stat_Pearson_Correlation_DDO1.txt",header = T, check.names = F)
    #MShannon_entropy_tab=read.table("./Table_Stat_Shannon_Entropy_DDO1.txt",header = T, check.names = F)
    #MRenyi_entropy_tab=read.table("./Table_Stat_MRenyi_Entropy_DDO1.txt",header = T, check.names = F)
    
    Mtaxa_tab$Tissue=gsub("_.*","",Mtaxa_tab$Samples)
    Mtaxa_tab = Mtaxa_tab[complete.cases(Mtaxa_tab), ]
    write.table(x = Mtaxa_tab, file = paste("./",destination_folder,"/","Table_Stat_Taxa_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    #Mtaxa_tab_mean =aggregate(Mtaxa_tab[,c(4)], list(Mtaxa_tab$Tissue,Mtaxa_tab$Samples,Mtaxa_tab$n), mean)
    #colnames(Mtaxa_tab_mean)=c("Tissue","Samples","n","Freq")
    #Mtaxa_tab_sd =aggregate(Mtaxa_tab[,c(4)], list(Mtaxa_tab$Tissue,Mtaxa_tab$Samples,Mtaxa_tab$n), sd)
    #colnames(Mtaxa_tab_sd)=c("Tissue","Samples","n","sd")
    #Mtaxa_tab1=merge(Mtaxa_tab_mean,Mtaxa_tab_sd,by=c("Tissue","Samples","n"),all = TRUE)
    
    Mcor=Mcor[!is.na(Mcor$id_pos),]
    write.table(x = Mcor, file = paste("./",destination_folder,"/","Table_Stat_Pearson_Correlation_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    MShannon_entropy_tab$Tissue=gsub("_.*","",MShannon_entropy_tab$pop)
    MShannon_entropy_tab = MShannon_entropy_tab[complete.cases(MShannon_entropy_tab), ]
    write.table(x = MShannon_entropy_tab, file = paste("./",destination_folder,"/","Table_Stat_Shannon_Entropy_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    MRenyi_entropy_tab$Tissue=gsub("_.*","",MRenyi_entropy_tab$pop)
    MRenyi_entropy_tab = MRenyi_entropy_tab[complete.cases(MRenyi_entropy_tab), ]
    write.table(x = MRenyi_entropy_tab, file = paste("./",destination_folder,"/","Table_Stat_MRenyi_Entropy_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    ner2 = ner2[complete.cases(ner2), ]
    #ner2=ner2[!duplicated( ner2[c("Samples")]),]
    ner2$se=ner2$sd/sqrt(length(CG_pos))
    ner2$Groups=sub("_[^_]+$", "", ner2$Samples)
    ner2$Groups=gsub("Random_control", "Random",  ner2$Groups)
    write.table(x = ner2, file = paste("./",destination_folder,"/","Table_Stat_Methylation_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
    
    ner2$Tissue=gsub("_.*","",ner2$Samples)
    
    Tissue= unique(ner2$Tissue)
    Tissue=gsub("Random$", '', Tissue)
    Tissue=stri_remove_empty(Tissue)
    
    # if option is TRUE Make a plot with All samples togheter ("TRUE" or "FALSE")
    if (Plot_All=="TRUE") {
      #Tissue=append(Tissue,"All")
      Tissue="All"
    } else {
      Tissue=Tissue
    }
    
    for (t in 1:length(Tissue)) {
      
      if (Tissue[t]=="All") {
        SamplesID=subset(map, map$Tissue!="Random")
        gdata=subset(ner2, ner2$Tissue!="Random")
      } else {
        SamplesID=map[map$Group %in% Tissue[t],]
        gdata=ner2[ner2$Tissue %in% Tissue[t],]
      }
      SamplesID=rbind(SamplesID,map[(dim(map)[1]),])
      SamplesID$order=paste( SamplesID$Group, SamplesID$Tissue,SamplesID$Description, SamplesID$Rep,sep="_")
      SamplesID$order[dim(SamplesID)[1]]="Random" #paste("Random_control_nCG",dim(CG_pos)[2],sep="")
      SamplesID$order1=gsub("^.*?\\_","", SamplesID$order)
      
      Random=subset(ner2, ner2$Tissue=="Random")
      #Remove all before first "_"
      Random$Samples=sub("_.*","",Random$Samples)
      Random$Timing=sub("_.*","",Random$Samples)
      
      #gdata=ner2[ner2$Tissue %in% Tissue[t],]
      
      gdata=gdata[match(SamplesID$order, gdata$Samples),]
      #Remove all after last"_"
      gdata$Timing=sub("[_][^_]+$","",gdata$Samples)
      gdata$Timing=sub(".*_","",gdata$Timing)
      
      gdata1=rbind(gdata,Random)
      gdata1=gdata1[!is.na( gdata1$Timing),]
      ########## MethylationPlotter for Loolipop
      tab=gdata1
      tab1=tab[,-c((length(CG_pos)+1):length(tab))]
      tab2=tab1
      tab2=format(round(tab2, 2), nsmall = 2)
      row.names(tab2)=tab$Samples
      #row.names(tab2)=paste(tab$Tissue,tab$Samples,sep="")
      tab2[dim(tab2)[1]+1,]=CG_pos[1,]
      row.names(tab2)[dim(tab2)[1]]="position"
      
      #Save and load on "http://maplab.imppc.org/methylation_plotter/"
      write.table(x=tab2, file=paste("./",destination_folder,"/","Methylation_Plotter_",Tissue[t],"_",gene, ".txt", sep=""), row.names=T, col.names = F,quote = F, sep="\t")
      
      ####### Fragment Methylation Plot
      #gdata1$Group=paste(tab$Tissue,tab$Samples,sep="_")
      gdata1$Group=gsub("^.*?\\_","", gdata1$Samples)
      if (Tissue[t]=="All") {
        gdata1$Group <- factor(gdata1$Group, levels = unique(gdata1$Group))
        gdata1$Timing <- factor(gdata1$Timing, levels = unique(gdata1$Timing))
        X=gsub("^.*?_","",gdata1$Group)
      } else {
        gdata1$Group <- factor(gdata1$Group, levels =  unique(gdata1$Group))
        gdata1$Timing <- factor(gdata1$Timing, levels = unique(gdata1$Timing))
        X=gsub("^.*?_","",gdata1$Group)
      }
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Meth_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Meth_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      j=ggplot(gdata1, aes(x=X, y=Means,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        geom_errorbar(aes(ymin=Means-se, ymax=Means+se)) +
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("%_Meth") +
        ggtitle(paste("Meth_",Tissue[t],"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      j=addSmallLegend(j)
      plot(j)
      dev.off()
      
      ####### Average Fragment Methylation Plot
      aggregateGroups=aggregate(gdata1[,(dim(CG_pos)[2]+1)], list(gdata1$Groups,gdata1$Tissue,gdata1$Timing), mean)
      colnames(aggregateGroups)=c("Samples","Tissue","Timing","mean")
      aggregateGroups_sd=aggregate(gdata1[,(dim(CG_pos)[2]+1)], list(gdata1$Groups,gdata1$Tissue,gdata1$Timing), sd)
      colnames(aggregateGroups_sd)=c("Samples","Tissue","Timing","sd")
      aggregateGroups_sd[is.na(aggregateGroups_sd)] <- 0
      aggregateGroups=merge(aggregateGroups,aggregateGroups_sd,by=c("Samples","Tissue","Timing"),all=TRUE)
      #Save Average Methylation
      write.table(x=aggregateGroups, file=paste("./",destination_folder,"/","Table_Stat_Average_Methylation_",Tissue[t],"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      aggregateGroups$Timing <- factor(aggregateGroups$Timing, levels = unique(map$Description))
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Average_Meth_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Average_Meth_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      jj=ggplot(aggregateGroups, aes(x=Timing, y=mean,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("%_Meth") +
        ggtitle(paste("Average_Meth_",Tissue[t],"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      jj=addSmallLegend(jj)
      plot(jj)
      dev.off()
      
      # Average CpGs Methylation Plot
      #ner3=ner2[ner2$Groups %in% gdata1$Groups,]
      ner3=melt(gdata1[,c((1:dim(CG_pos)[2]),(dim(gdata1)[2]-5),(dim(gdata1)[2]-2))], id=c("Samples","Tissue"))
      colnames(ner3)=c("Samples","Tissue","CpGs","Freq")
      #ner3$Rep=sub(".*_","",ner3$Samples)
      ner3$Samples=sub("_M.*","",ner3$Samples)
      ner3$Timing=sub(".*_","",ner3$Samples)
      #ner3$Samples=paste(ner3$Tissue,ner3$Samples,sep="_")
      
      ner3$Timing <- factor(ner3$Timing, levels = unique(map$Description))
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","CpGs_Meth_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","CpGs_Meth_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      p<-ggplot(ner3, aes(x=CpGs, y=Freq, group=Timing)) +
        geom_line(aes(color=Timing))+
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        ggtitle(paste("CpGs_Meth_",Tissue[t],"_",gene,sep="")) + 
        geom_point(aes(color=Timing))+
        theme_classic() +
        scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      # theme(legend.position="bottom")
      
      # Apply on original plot
      p=addSmallLegend(p)
      #p+guides(fill=guide_legend(ncol=2))
      #p+ guides(col = guide_legend(nrow = unique(dim(ner3$Samples))))
      plot(p)
      dev.off()
      
      ####### Average CpGs Methylation Plot
      aggregateCpGsGroups=aggregate(ner3[,4], list(ner3$Samples,ner3$Tissue,ner3$CpGs), mean)
      colnames(aggregateCpGsGroups)=c("Samples","Tissue","CpGs","mean")
      aggregateCpGsGroups_sd=aggregate(ner3[,4], list(ner3$Samples,ner3$Tissue,ner3$CpGs), sd)
      colnames(aggregateCpGsGroups_sd)=c("Samples","Tissue","CpGs","sd")
      
      std <- function(x) sd(x)/sqrt(length(x))
      aggregateCpGsGroups_se=aggregate(ner3[,4], list(ner3$Samples,ner3$Tissue,ner3$CpGs), std)
      colnames(aggregateCpGsGroups_se)=c("Samples","Tissue","CpGs","se")
      
      aggregateCpGsGroups=merge(aggregateCpGsGroups,aggregateCpGsGroups_sd,by=c("Samples","Tissue","CpGs"),all=TRUE)
      aggregateCpGsGroups=merge(aggregateCpGsGroups,aggregateCpGsGroups_se,by=c("Samples","Tissue","CpGs"),all=TRUE)
      aggregateCpGsGroups$Timing=sub(".*_","",aggregateCpGsGroups$Samples)
      #Save Average Methylation
      write.table(x=aggregateCpGsGroups, file=paste("./",destination_folder,"/","Table_Stat_Average_Methylation_",Tissue[t],"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Average_CpGs_Meth+se_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Average_CpGs_Meth+se_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      p1<-ggplot(aggregateCpGsGroups, aes(x=CpGs, y=mean,group=Samples,color=Timing)) +
        geom_line()+
        geom_point()+
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
        facet_wrap(~ Tissue, scales="free_x") +
        ylab("Average_CpGs_Meth+se")+
        ylim(0,1)+
        ggtitle(paste("Average_CpGs_Meth+se_",Tissue[t],"_",gene,sep="")) + 
        geom_point(aes(color=Timing))+
        theme_classic() +
        scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
      
      # Apply on original plot
      p1=addSmallLegend(p1)
      plot(p1)
      dev.off()
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Average_CpGs_Meth_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Average_CpGs_Meth_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      p2<-ggplot(aggregateCpGsGroups, aes(x=CpGs, y=mean,group=Samples,color=Timing)) +
        geom_line()+
        geom_point()+
        facet_wrap(~ Tissue, scales="free_x") +
        ylab("Average_CpGs_Meth")+
        ylim(0,1)+
        ggtitle(paste("Average_CpGs_Meth_",Tissue[t],"_",gene,sep="")) + 
        geom_point(aes(color=Timing))+
        theme_classic() +
        scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
      
      # Apply on original plot
      p2=addSmallLegend(p2)
      plot(p2)
      dev.off()
      
      rm(ner3,gdata1,aggregateCpGsGroups,aggregateGroups,Random)
      
      ############################## Taxa Distribution 
      Random=subset(Mtaxa_tab, Mtaxa_tab$Tissue=="Random")
      
      if (Tissue[t]=="All") {
        taxa=subset(Mtaxa_tab, Mtaxa_tab$Tissue!="Random")
        #taxa=Mtaxa_tab
      } else {
        taxa=Mtaxa_tab[Mtaxa_tab$Tissue %in% Tissue[t],]
      }
      
      taxa1=rbind(taxa,Random)
      
      class="Methylated"
      
      #Remove NA
      taxa1=taxa1[!is.na(taxa1$Samples),]
      #taxa1 = taxa[complete.cases(taxa), ]
      write.table(x = taxa1, file = paste("./",destination_folder,"/","Table_Taxonomy_",Tissue[t],"_",class,"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      #taxa=read.table(paste("./Table_Taxonomy_",Tissue[t],"_",class,"_",gene, ".txt",sep=""),header = T, check.names = F)
      #taxa=taxa[,map$id] 
      taxa1$n=factor(taxa1$n)
      taxa1$Samples=gsub("^.*?\\_","", taxa1$Samples)
      taxa1$Samples=factor(taxa1$Samples,levels=unique(taxa1$Samples))
      taxa1$Samples=gsub("^.*?_","",taxa1$Samples)
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Taxonomy_",Tissue[t],"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Taxonomy_",Tissue[t],"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      p3<-ggplot(taxa1, aes(x=Samples, y=Freq,fill = forcats::fct_rev(n))) +
        geom_bar(stat="identity", color="black",position="fill")+
        geom_text(aes(label=round(Freq,digits = 2)),size = 2,color = "black",position=position_fill(vjust=0.5))+
        #geom_text(aes(label = ),  size = 3, position = position_dodge(width = 0.9),vjust = -2)
        facet_wrap(~ Tissue, scales="free_x") +
        #ylab("Freq of Species")+
        xlab("Samples")+
        #ylim(0,1)+
        ggtitle(paste("Taxonomy_",Tissue[t],"_",gene,sep="")) + 
        theme_classic() +
        labs(fill="Species")+
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        guides(fill=guide_legend(ncol=1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "right", legend.box = "vertical")
      
      # Apply on original plot
      p3=addSmallLegend(p3)
      plot(p3)
      dev.off()
      
      rm(taxa1,Random)
      ################## Samples Correlation
      cor1=Mcor[,SamplesID$order]
      cor1=cbind(cor1,Mcor$id_pos)
      colnames(cor1)[(dim(cor1)[2])]="id_pos"
      cor1=cor1[,c("id_pos",SamplesID$order)]
      #cor=cbind(cor1,Random)
      cor=cor1
      class="Methylated"
      
      cor[is.na(cor)] <- 0
      
      cor=data.matrix(cor[,c(2:dim(cor)[2]) ])
      
      cor2<-rcorr(cor)
      cor3=flattenCorrMatrix(cor2$r, cor2$P)
      write.table(x=cor3, file=paste("./",destination_folder,"/","Table_Pearson_Correlation_matrix_",Tissue[t],"_",class,"_",gene,".txt", sep=""), row.names=T, col.names = F,quote = F, sep="\t")
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Pearson_Correlation_matrix_",Tissue[t],"_",class,"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Pearson_Correlation_matrix_",Tissue[t],"_",class,"_",gene,".png", sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }

      #Fix sample name
      colnames(cor)=gsub(paste("^.*?",Tissue[t],"_",sep=""),"",colnames(cor))
      
      chart.Correlation(cor,histogram = TRUE,method ="pearson")
      mtext(paste("Pearson_Correlation_matrix_",Tissue[t],"_",class,"_",gene, sep=""), outer=TRUE,  cex=1, line=-1)
      #plot(chart1)
      dev.off()
      
      ########################## Shannon Entropy Plot 
      
      if (Tissue[t]=="All") {
        Shannon_Entropy=subset(MShannon_entropy_tab,MShannon_entropy_tab$Tissue!="Random")
      } else {
        Shannon_Entropy= MShannon_entropy_tab[MShannon_entropy_tab$Tissue %in% Tissue[t],]
      }
      
      write.table(x = MShannon_entropy_tab, file = paste("./",destination_folder,"/","Table_Stat_Shannon_Entropy_",Tissue[t],"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      Random=subset(MShannon_entropy_tab, MShannon_entropy_tab$Tissue=="Random")
      #Remove all before first "_"
      Random$pop=sub("_.*","",Random$pop)
      Random$Timing=sub("_.*","",Random$pop)
      
      Shannon_Entropy=Shannon_Entropy[match(SamplesID$order, Shannon_Entropy$pop),]
      #Remove all after last"_"
      Shannon_Entropy$Timing=sub("[_][^_]+$","",Shannon_Entropy$pop)
      Shannon_Entropy$Timing=sub(".*_","",Shannon_Entropy$Timing)
      
      Shannon_Entropy1=rbind(Shannon_Entropy,Random)
      #Remove NA
      Shannon_Entropy1=Shannon_Entropy1[!is.na(Shannon_Entropy1$pop),]
      #Shannon_Entropy1=Shannon_Entropy
      
      
      Shannon_Entropy1$Group=gsub("^.*?\\_","",  Shannon_Entropy1$pop)
      #Shannon_Entropy1$Group=paste(Shannon_Entropy1$Tissue,Shannon_Entropy1$pop,sep="_")
      if (Tissue[t]=="All") {
        #Shannon_Entropy1$Group <- factor(Shannon_Entropy1$Group, levels = Shannon_Entropy1$Group)
        #X=Shannon_Entropy1$Group
        #Shannon_Entropy1$pop <- factor(Shannon_Entropy1$pop, levels = Shannon_Entropy1$pop)
        Shannon_Entropy1$Group <- factor(Shannon_Entropy1$Group, levels = unique(Shannon_Entropy1$Group))
        Shannon_Entropy1$Timing <- factor(Shannon_Entropy1$Timing, levels = unique(map$Description))
        X=Shannon_Entropy1$Timing
      } else {
        #Shannon_Entropy1$pop <- factor(Shannon_Entropy1$pop, levels = Shannon_Entropy1$pop)
        Shannon_Entropy1$Group <- factor(Shannon_Entropy1$Group, levels = unique(Shannon_Entropy1$Group))
        Shannon_Entropy1$Timing <- factor(Shannon_Entropy1$Timing, levels = unique(map$Description))
        X=Shannon_Entropy1$Timing
      }
      
      
      class="Methylated"
      
      write.table(x = Shannon_Entropy1, file = paste("./",destination_folder,"/","Table_Shannon_Entropy_",Tissue[t],"_",class,"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Shannon_Entropy_",Tissue[t],"_",class,"_",gene,"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Shannon_Entropy_",Tissue[t],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
      }
      j1=ggplot(Shannon_Entropy1, aes(x=X, y=shannon_entropy1,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("Shannon_Entropy") +
        ggtitle(paste("Shannon_Entropy_", Tissue[t],"_",class,"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        #guides(fill=guide_legend(ncol=1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      j1=addSmallLegend(j1)
      plot(j1)
      dev.off()
      
      #Average Shannon
      aggregateGroups=aggregate(Shannon_Entropy1[,2], list(Shannon_Entropy1$Timing,Shannon_Entropy1$Tissue), mean)
      colnames(aggregateGroups)=c("Timing","Tissue","mean")
      aggregateGroups_sd=aggregate(Shannon_Entropy1[,2], list(Shannon_Entropy1$Timing,Shannon_Entropy1$Tissue), sd)
      colnames(aggregateGroups_sd)=c("Timing","Tissue","sd")
      aggregateGroups=merge(aggregateGroups,aggregateGroups_sd,by=c("Timing","Tissue"),all=TRUE)
      #Save Average Methylation
      write.table(x=aggregateGroups, file=paste("./",destination_folder,"/","Table_Stat_Average_Shannon_Entropy_",Tissue[t],"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      aggregateGroups$Timing <- factor(aggregateGroups$Timing, levels = unique(map$Description))
      
      #Average Shannon plot
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Average_Shannon_Entropy_", Tissue[t],"_",class,"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Average_Shannon_Entropy_", Tissue[t],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      jj=ggplot(aggregateGroups, aes(x=Timing, y=mean,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("Average_Shannon_Entropy") +
        ggtitle(paste("Average_Shannon_Entropy_", Tissue[t],"_",class,"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        #guides(fill=guide_legend(ncol=1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      jj=addSmallLegend(jj)
      plot(jj)
      dev.off()
      
      rm(aggregateGroups,Shannon_Entropy1,Random)
      
      ########################## Renyi Entropy Plot
      
      Random=subset(MRenyi_entropy_tab, MRenyi_entropy_tab$Tissue=="Random")
      #Remove all before first "_"
      Random$pop=sub("_.*","",Random$pop)
      Random$Timing=sub("_.*","",Random$pop)
      
      if (Tissue[t]=="All") {
        Renyi_Entropy=subset(MRenyi_entropy_tab,MRenyi_entropy_tab$Tissue!="Random")
      } else {
        Renyi_Entropy=MRenyi_entropy_tab[MRenyi_entropy_tab$Tissue %in% Tissue[t],]
      }
      
      Renyi_Entropy=Renyi_Entropy[match(SamplesID$order,  Renyi_Entropy$pop),]
      #Remove all after last"_"
      Renyi_Entropy$Timing=sub("[_][^_]+$","",Renyi_Entropy$pop)
      Renyi_Entropy$Timing=sub(".*_","",Renyi_Entropy$Timing)
      
      Renyi_Entropy1=rbind(Renyi_Entropy,Random)
      #Renyi_Entropy1=Renyi_Entropy
      #Remove NA
      Renyi_Entropy1=Renyi_Entropy1[!is.na(Renyi_Entropy1$pop),]
      
      Renyi_Entropy1$Group=gsub("^.*?\\_","",  Renyi_Entropy1$pop)
      #Renyi_Entropy1$Group=paste(Renyi_Entropy1$Tissue,Renyi_Entropy1$pop,sep="_")
      if (Tissue[t]=="All") {
        #Renyi_Entropy1$Group <- factor(Renyi_Entropy1$Group, levels = Renyi_Entropy1$Group)
        #X=Renyi_Entropy1$Group
        Renyi_Entropy1$Group <- factor(Renyi_Entropy1$Group, levels = unique(Renyi_Entropy1$Group))
        Renyi_Entropy1$Timing <- factor(Renyi_Entropy1$Timing, levels = unique(map$Description))
        X=Renyi_Entropy1$Timing
      } else {
        Renyi_Entropy1$Group <- factor(Renyi_Entropy1$Group, levels = unique(Renyi_Entropy1$Group))
        Renyi_Entropy1$Timing <- factor(Renyi_Entropy1$Timing, levels = unique(map$Description))
        X=Renyi_Entropy1$Timing
      }
      
      class="Methylated"
      
      write.table(x = Renyi_Entropy1, file = paste("./",destination_folder,"/","Table_Renyi_Entropy_",Tissue[t],"_",class,"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Renyi_Entropy_",Tissue[t],"_",class,"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Renyi_Entropy_",Tissue[t],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
      }
      j2=ggplot(Renyi_Entropy1, aes(x=X, y=RenyiEntropy,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("Renyi_Entropy") +
        ggtitle(paste("Renyi_Entropy_", Tissue[t],"_",class,"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        #guides(fill=guide_legend(ncol=1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      j2=addSmallLegend(j2)
      plot(j2)
      dev.off()
      
      #Average Renyi
      aggregateGroups=aggregate(Renyi_Entropy1[,2], list(Renyi_Entropy1$Timing,Renyi_Entropy1$Tissue), mean)
      colnames(aggregateGroups)=c("Timing","Tissue","mean")
      aggregateGroups_sd=aggregate(Renyi_Entropy1[,2], list(Renyi_Entropy1$Timing,Renyi_Entropy1$Tissue), sd)
      colnames(aggregateGroups_sd)=c("Timing","Tissue","sd")
      aggregateGroups=merge(aggregateGroups,aggregateGroups_sd,by=c("Timing","Tissue"),all=TRUE)
      #Save Average Methylation
      write.table(x=aggregateGroups, file=paste("./",destination_folder,"/","Table_Stat_Average_Renyi_Entropy_",Tissue[t],"_",gene, ".txt",sep=""), quote=F,sep="\t",row.names=F)
      
      aggregateGroups$Timing <- factor(aggregateGroups$Timing, levels = unique(map$Description))
      
      #Average Shanno plot
      if (saveAs=="pdf") {
        pdf(paste("./",destination_folder,"/","Average_Renyi_Entropy_", Tissue[t],"_",class,"_",gene,".pdf", sep=""),paper = "a4")
      } else {
        png(paste("./",destination_folder,"/","Average_Renyi_Entropy_", Tissue[t],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8)
      }
      jj=ggplot(aggregateGroups, aes(x=Timing, y=mean,fill=Timing)) + 
        geom_bar(position=position_dodge(), stat="identity",
                 colour="black", # Use black outlines,
                 size=.3) +      # Thinner lines
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
        facet_wrap(~ Tissue, scales="free_x") +
        ylim(0,1)+
        xlab("Samples") +
        ylab("Average_Renyi_Entropy") +
        ggtitle(paste("Average_Renyi_Entropy_", Tissue[t],"_",class,"_",gene,sep="")) + 
        theme_classic() +
        #scale_fill_continuous(guide = guide_legend()) +
        theme(legend.position="right", legend.direction="vertical")+
        #guides(fill=guide_legend(ncol=1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"))
      
      addSmallLegend <- function(myPlot, pointSize = 6, textSize = 6, spaceLegend = 0.1) {
        myPlot +
          guides(shape = guide_legend(override.aes = list(size = pointSize)),
                 color = guide_legend(override.aes = list(size = pointSize))) +
          theme(legend.title = element_text(size = textSize), 
                legend.text  = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
      }
      
      # Apply on original plot
      jj=addSmallLegend(jj)
      plot(jj)
      dev.off()
      
      rm(aggregateGroups,Renyi_Entropy1,Random)
    }
    
    
    #Delete everything and leave only the.out files
    # rm(list = ls()[!ls() %in% c("saveAs","Meth_unMeth_corr","dectobin","chi_square","flattenCorrMatrix","radian.rescale","quantile","nProcessor","minimal_common_structure","cores_status_assignment",
    #                             "cores_heatmap","epialleles_heatmap","strings_to_BinaryProfiles","complexity","compare_Samples_epialleles","combs_summ_clust1",
    #                             "combs_summary1","combinations","summ_stat","make_graph","cores_extraction","tab_meth","Mtab_clones","co_Meth_occurence","tetrachoric_corr",
    #                             "nProcessor","genes","maps","Complexity","y","pvalue","Reads","Remove","min","maxCG_pos","class", 
    #                             "co_Meth_unMeth_occurence","co_unMeth_Meth_occurence","correlation_plot","Plot_All")])
    
    rm(list = ls()[!ls() %in% c("splitForGene","destination_folder","Modify_mono","saveAs","Meth_unMeth_corr","Plot_All","nProcessor","genes","maps","y","pvalue","Reads","Remove","min","maxCG_pos","class")])
}

