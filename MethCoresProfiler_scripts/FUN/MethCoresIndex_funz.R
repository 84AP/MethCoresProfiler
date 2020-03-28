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

MethCoresIndex <- function(timing,t,tdim5,tmap1,g,timing1,t1) {
  
    my_data = subset(tdim5, Group ==timing[t])
    my_data$Rep=gsub(".*_","",my_data$Samples) 
    samples=as.vector(unique(my_data$Group))
    
    if (dim(my_data)[1]==0) {
      my_data2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
      cores1=""
    } else {
    
    for (i in 1:length(cores_cutOff)) {
      #Use minimal_common_structure2==not take common epiallels amon replicate but only common CpGs among all epiallels
    cores=minimal_common_structure2(cor=my_data,tmap1=tmap1,timing=timing,t=t,v=i)
    cores1=unlist(strsplit(cores$cores, "-"))
    
    if (length(cores1)<=1) next
    my_data$quantile=cores$quantile
    #save list epialleles group
    write.table(x=my_data, file=paste("./",destination_folder,"/Table_Significant_Epialleles_",Groups[g],"_",timing[t],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T,col.names=T, quote = F, sep="\t")
    
    my_data1=aggregate( Freq ~ Rep,  my_data, sum )
    my_data1$Group=unique(my_data$Group)
    my_data1$core=paste(cores1,collapse = "-")
    
    #my_data2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
    #my_data2=my_data2[0,]
   Freq_core=mean(my_data1$Freq)
    
    #Crea lista di tutte le CpG negli epialleli significativi
    cpg_sets1 =unlist(strsplit(my_data$core,"-"))
    tab_cpg_sets1=as.data.frame(table(cpg_sets1))
    tab_cpg_sets1=tab_cpg_sets1[!tab_cpg_sets1$cpg_sets1 %in% cores1, ]
    tab_cpg_sets1$cpg_sets1=as.character(tab_cpg_sets1$cpg_sets1)
    tab_cpg_sets1$Freq=as.numeric(tab_cpg_sets1$Freq)
    sum_Freq_cpg_set1=as.numeric(sum(tab_cpg_sets1$Freq))
    
    #for (mm in 1:nrow(my_data1)) {
    my_data2=data.frame(matrix(NA, nrow = length(setdiff(unlist(CG_pos[1,]),cores1)), ncol = 5))
    colnames(my_data2)=c("From","to","weight","quantile","sample")
    my_data2$From=unique(my_data1$core)
    my_data2$quantile=cores$quantile
    my_data2$sample=samples
    
    if (dim(tab_cpg_sets1)[1]==0) {
      my_data2$to[1]= ""
      my_data2$weight[1]=  Freq_core 
    } else {
    for (n in 1:nrow(tab_cpg_sets1)) {
      tab_cpg_sets1[n,"Freq"]=tab_cpg_sets1[n,"Freq"]/ sum_Freq_cpg_set1
      my_data2$to[n]= tab_cpg_sets1$cpg_sets1[n]
    my_data2$weight[n]= Freq_core*tab_cpg_sets1[n,"Freq"]
    }
      }   
    my_data2=na.omit(my_data2)
    #mydata2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
      #mydata1=my_data1[mm,]
      
      #crea una lista con il set di cpg per ogni epiallele (la puoi mettere pure fuori al loop questa)
      #cpg_sets1 =unlist(strsplit(mydata1$core,"-"))
      #cpg_sets1 =unlist(strsplit(my_data$core,"-"))
      
     
      
    #   #cores1
    #   #if  (all(cores1 %in% cpg_sets1)==TRUE) {
    #     mydata2$From[1]=paste(cores1,collapse = "-")
    #     mydata2$to[1]=paste(setdiff(cpg_sets1,cores1),collapse = "-")
    #     mydata2$weight[1]=mydata1$Freq
    #     mydata2$quantile[1]=cores$quantile
    #     mydata2$sample=samples
    #     #mydata2=mydata2[!(is.na(mydata2$From) | mydata2$From==""), ]
    #   }  else {
    #     mydata2$From=paste(cpg_sets1,collapse = "-")
    #     mydata2$to=""
    #     mydata2$weight=mydata1$Freq
    #     mydata2$quantile[1]=""
    #     mydata2$sample=samples
    #   }
    #   #replace blank with cores
    #   my_data2=rbind(my_data2,mydata2)
    #  
    # }
    # 
    # my_data2=my_data2[!(is.na(my_data2$From) | my_data2$From==""), ]
    
    if (Groups[g]=="All") {
      my_data2$Tissue=timing1[t1]
    } else {
      my_data2$Tissue=timing[t]
    }
    
    #Change cores dimension
    #if (sum(my_data2$to != "")<1) next
    #if (sum(my_data2$to != "")>=1) break;

    if (sum(my_data2$From != "")>=1) break;

    }

    my_data2=list(my_data=my_data,my_data1=my_data1,my_data2=my_data2,cores=cores1)
    
    }
  return(my_data2)
}
