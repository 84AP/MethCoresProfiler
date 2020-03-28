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

MethCores_tab2 <- function(timing,t,tdim5,tmap1,g,filename,Complexity,hh,timing1,t1,Epi_list1) {

    MethCores_Index=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Group=NA,Groups=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1)
  tab_Stat_MethCores=data.frame(Samples=NA,average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,MethCores=NA,average_MethCores_Index=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,average_Clonality_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1)
  tab_my_data= data.frame(Rep=NA,Freq=-1,Group=NA,core=NA)
  
  my_data2=MethCoresIndex(timing=timing,t=t,tdim5,tmap1=tmap1,g,timing1,t1)
  #replace blank with cores
  # my_data2$From <- sub("^$", cores[1],  my_data2$From)
  cores=my_data2$cores
  my_data3=my_data2$my_data2
  #Riempi caselle vuote con colonno adiacenti
  my_data3[] <- lapply(my_data3, as.character)
  
  my_data1=my_data2$my_data1
  
  #my_data2$to[is.na(my_data2$to)] <- as.character(my_data2$From[is.na(my_data2$to)])
  
 # if (unique(grepl('_',tmap1$SamplesID)[-length(unique(tmap1$SamplesID))])) {
  
  # if (Groups[g]=="All") {
  #   my_data3$Group=paste(timing1[t1],timing[t],sep="_")
  # } else {
  #   my_data3$Group=paste(Groups[g],timing[t],sep="_")
  # }
  # } else {
    my_data3$Group=timing[t]
  #}
  
  my_data3$Complexity=Complexity[hh]
  my_data3$Tissue=Groups[g]
  
  if (population_weight==TRUE) {
    #Load Population
    population=Epi_list1[grep(timing[t], Epi_list1, fixed=T)]
    population1=read.table(population, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    cpg_sets=strsplit(population1$id_pos,"-")
    core_index=which(sapply(cpg_sets,function(x) { all(cores %in% x )}))
    population2= population1[core_index,]
    my_data3$weight=sum(population2$FreqCloneForSample)
  } else {
    my_data3$weight <- as.numeric(as.character(my_data3$weight))
  }
  
  #if (sum(is.na(my_data3$From))==dim(my_data3)[1]) next
  
  if (sum(my_data3$to != "")>1) {
    my_data3$MethCores <- paste(cores,collapse = "-")
    #my_data3=my_data3[!(is.na(my_data3$to) | my_data3$to==""), ]
  } else {
    my_data3=my_data3[order(-my_data3$weight),]
    my_data3$MethCores <- my_data3$From[1]
  }
  
  #Remove empty row
  #my_data3=my_data3[!(is.na(my_data3$From) | my_data3$From==""), ]
  
  #Assign colours amd size to each vertex
  From=subset(my_data3, select=-c(to))
  colnames(From)[1]="core"
  to=subset(my_data3, select=-c(From))
  colnames(to)[1]="core"
  df=rbind(From,to)
  df$weight <- as.numeric(as.character(df$weight))
  #Remove empty row
  df=df[!(is.na(df$core) | df$core==""), ]
  df$colors=tdim5$Selection[match(df$MethCores, tdim5$core)]
  ##Remove duplicated row
  df=df %>% distinct(core, weight, .keep_all = TRUE)
  vertex=aggregate(df$weight, by=list(Category=df$core),FUN=sum)
  colnames(vertex)=c("vertex","Freq")
  rownames(vertex)=1:dim(vertex)[1]
  #vertex$colors=tdim5$Selection[match(vertex$vertex, tdim5$core)]
  #vertex$colors[is.na(vertex$colors)] <- "blue"
  vertex$colors= ifelse(vertex$vertex == (paste(cores,collapse = "-")), "red", "blue")
  
  # if there are is not red take a more frequently
  #if ("red" %in% vertex$colors) {
  #  vertex$Selection=vertex$colors
  #} else {
   # vertex$Selection= ifelse(vertex$vertex == unique(my_data3$MethCores), "red", "blue")
  #}
  
  #Remove all less to 0.1%
  rvertex=vertex[!(vertex$vertex==""),]
  if (max(rvertex$Freq)>=vertex_cutOFF) {
    #my_data3=my_data3[!my_data3$weight <= vertex_cutOFF,]
    vertex=vertex[!vertex$Freq <= vertex_cutOFF,]
  } else {
    #my_data3=my_data3[order(-my_data3$weight),]
    #my_data3 <- my_data3[1,]
    #vertex=vertex[!is.na(vertex$vertex),]
    vertex=vertex[order(-vertex$Freq),]
    vertex <- vertex[1:2,]
  }
  vertex=vertex[!is.na(vertex$vertex),]
  
  #Replace empty with "NA"
  #my_data3[my_data3==""] <- "NA"
  rownames(my_data3)=1:dim(my_data3)[1]
  #vertex1=my_data3[!grepl("", my_data3$to),]
  if (max(my_data3$weight)>=vertex_cutOFF) {
    #vertex1=my_data3[!(is.na(my_data3$to) | my_data3$to==""), ]
    vertex1=data.frame(From=vertex$vertex,to="",weight=0)
  } else {
    vertex1=data.frame(From=vertex$vertex,to="",weight=0)
  }
  #vertex1=vertex1[!(vertex1$From==""),]
  vertex1[] <- lapply(vertex1, as.character)
  vertex1$weight=as.numeric(vertex1$weight)
  
  tab_vertex1=setdiff(vertex$vertex,vertex1$From)
  tab_vertex1=stri_remove_empty_na(tab_vertex1)
  #remove vertex from vertex1 not present in vertex
  vertex1=vertex1[vertex1$to %in% tab_vertex1,]
  
  #vertex1=my_data3
  write.table(x= my_data3, file=paste("./",destination_folder,"/Table_Network_Epialleles_",filename,"_",timing[t],"_",class,"_",gene, ".txt", sep=""), row.names=F,col.names=T, quote = F, sep="\t")
  
  # create data:
  links=vertex1[,1:2]
  
  # Turn it into igraph object
  network=graph_from_data_frame(links, directed=F,vertices = vertex$vertex) 
  
  if (any(vertex$vertex=="")==TRUE) {
    network= delete_vertices(network, c("")) # here's my condition.
  } else {
    network=network
  }
  
  sz <- vertex$Freq[match(V(network)$name, vertex$vertex)]
  colors <- vertex$colors[match(V(network)$name, vertex$vertex)]
  
  #V(network)$name=vertex[,'vertex']
  V(network)$size=(sz*10)
  V(network)$color=colors
  E(network)$weight=(vertex1[,'weight']*10)
  #V(network)$edge=vertex1[,1:2]
  
  # Count the number of degree for each node:
  #deg=as.data.frame(degree(network, mode="all"))
  # Fixing ego
  minC <- rep(-Inf, vcount(network))
  maxC <- rep(Inf, vcount(network))
  minC[1] <- maxC[1] <- 0
  co <- layout_with_fr(network, minx=minC, maxx=maxC,
                       miny=minC, maxy=maxC)
  # Plot
  png(paste("./",destination_folder,"/Cores_graph_",filename,"_",timing[t],"_",class,"_",gene,".png",sep=""),
      width = 5*150,        # 5 x 300 pixels
      height = 5*150,
      res = 150,            # 300 pixels per inch
      pointsize = 4)        # smaller font size
  
  # add the size and border colour vectors to relevant plot arguments     
  plot(network, layout=co,vertex.size=V(network)$size,vertexedge=V(network)$edge,vertex.color=V(network)$colors, 
       edge.width=E(network)$weight,vertex.frame.color= colors, vertex.frame.width=2,vertex.label.dist=3,
       vertex.label.cex=1.5, add=FALSE,vertex.label.font=1,margin=0.5)
  title(paste("Cores_graph_",filename,"_",timing[t],"_",class,"_",gene,sep=""),cex.main=2,col.main="black") #
  
  #Legenda Interactions
  sizeCut1<- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  sizeCutScale1 <- sizeCut1
  
  a=legend('bottomleft',title = paste("Frequency of Interactions",sep=""),legend=unique(sizeCut1),lwd=sizeCutScale1,col='black',cex=1)
  #a <- legend('bottomleft',title = paste("Frequency of Interactions",sep=""),legend=unique(sizeCut),lwd=sizeCutScale/200,col='white', cex=1)
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  
  #legdata=legend('topleft',title = paste("Freq of Interactions",sep=""),legend=levels(sizeCut),pt.cex=(scaled/ns),col='black',lty=1:1, cex=1)
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Freq of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Freq of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  
  #vertex=vertex[!(is.na(vertex$vertex) | vertex$vertex==""), ]
  
  #Legend Cores size
  sizeCut<- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  sizeCutScale <- sizeCut*10
  
  legend('topleft',title = paste("Freq of Cores Occurence",sep=""),legend=unique(sizeCut),pt.cex= sizeCutScale,col='black')
  a <- legend('topleft',title = paste("Freq of Cores Occurence",sep=""),legend=unique(sizeCut),pt.cex=sizeCutScale/200,col='white', pch=21, pt.bg='white')
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  symbols(x,y,circles=sizeCutScale/200,inches=FALSE,add=TRUE,bg='orange')
  #legdata=legend('bottomright',title = paste("Freq of Cores Occurence",sep=""),legend=levels(sizeCut1),pt.cex=(scaled1/ns1),col='black',pch=21)
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Freq of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Freq of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  
  #Legenda Selection
  cores_cutOff1=unique(my_data3$quantile)
  cores_cutOff1= cores_cutOff1[ cores_cutOff1 != ""]
  legdata1=legend('bottomright',title = "Types of Selection",legend=c(paste("Common core >=",cores_cutOff1,"Percentile,pvalue<=",pvalue,sep=""),paste("UnCommon core,pvalue<=",pvalue,sep="")),col=c("red", "blue"), pch=21, cex=1.2)
  #legend(list(x=legdata1$rect$left,y=legdata1$rect$top+0.5),legend=c("Negative Selection, pvalue<=0.01","No Selection","Positive Selection, pvalue<=0.01"),col=c("blue","gray", "red"), lty=1:1, cex=0.8)
  #legend(list(x=legdata1$rect$left,y=legdata1$rect$top+0.5),legend=c("Negative Selection, pvalue<=0.01","No Selection","Positive Selection, pvalue<=0.01"),col=c("blue","gray", "red"), lty=1:1, cex=0.8)
  
  
  
  if (Groups[g]=="All") {
    samples=unique(my_data3$Group)
    samples= na.omit(samples)
    Stat$Group=sub(".*_","",Stat$Groups)
  } else {
    samples=unique(my_data3$Group)
    samples= na.omit(samples)
    Stat$Group=sub(".*_","",Stat$Groups)
    #Stat$id=paste(Stat$Tissue,Stat$Group,sep="_")
  }
  Stat$id=paste(Stat$Tissue,Stat$Group,sep="_")
  
  Stat1=Stat[Stat$Groups %in% samples,]
  # if there are is not red take a more frequently
  if ("red" %in% vertex$colors) {
    Stat_cores=vertex[vertex$colors == "red",]
  } else {
   vertex=vertex[order(-vertex$Freq),]
   Stat_cores= vertex[1,]
  }
  
  # Create list Samples
  Samples=as.vector(unique(Stat1$Samples))
  tab_Stat_cores=data.frame(matrix(ncol=3,nrow=length(Samples)))
  colnames( tab_Stat_cores)=c("Samples","average_methylation","MethCores_Index")
  
  my_data1=my_data2$my_data1
  
  for (dd in 1:length(Samples)) {
    ddim5=tdim5[paste(tdim5$Tissue,tdim5$Samples,sep="_") %in% Samples[dd], ]
    #ddim5=tdim5[paste(tdim5$Tissue,tdim5$Samples,sep="_") %in% Samples[dd], ]
    #check if core is present in string
    #grepl(Stat_cores$vertex, ddim5$core)
    ddim6= ddim5[grep(Stat_cores$vertex, ddim5$core), ]
    #ddim6= ddim5[ ddim5$core %in% Stat_cores$vertex, ]
    tab_Stat_cores$Samples[dd]=Samples[dd]
    tab_Stat_cores$MethCores_Index[dd]=my_data1$Freq[match(tab_Stat_cores$Samples[dd],paste(select,my_data1$Group,my_data1$Rep,sep="_"))]
    tab_Stat_cores$average_methylation[dd]=Stat1$Means[match(tab_Stat_cores$Samples[dd],Stat1$Samples)]
    #tab_Stat_cores$MoreRepresentativeEpiallele[dd]=sum(ddim6$Freq)
    
  }
  
  tab_Stat_cores$Clonality_Index=tab_Stat_cores$MethCores_Index/tab_Stat_cores$average_methylation
  tab_Stat_cores$MethCores=Stat_cores$vertex
  #Stat MethCores_Index
  
  stat_MethCores_Index <- stat.desc(tab_Stat_cores$MethCores_Index) 
  #stat_MethCores_Index=describe.by(tab_Stat_cores$MethCores_Index)
  # replace NA with 0 if there are not replicate
  #stat_MethCores_Index[is.na( stat_MethCores_Index)] <- 0
  
  tab_Stat_cores$average_MethCores_Index= stat_MethCores_Index["mean"]
  tab_Stat_cores$sd_MethCores_Index= stat_MethCores_Index["std.dev"]
  tab_Stat_cores$se_MethCores_Index= stat_MethCores_Index["SE.mean"]
  tab_Stat_cores=tab_Stat_cores[!is.na(tab_Stat_cores$Samples),]
  
  #Stat Clonality_Index
  stat_Clonality_Index=stat.desc(tab_Stat_cores$Clonality_Index)
  # replace NA with 0 if there are not replicate
  #stat_Clonality_Index[is.na(stat_Clonality_Index)] <- 0
  
  tab_Stat_cores$average_Clonality_Index= stat_Clonality_Index["mean"]
  tab_Stat_cores$sd_Clonality_Index= stat_Clonality_Index["std.dev"]
  tab_Stat_cores$se_Clonality_Index= stat_Clonality_Index["SE.mean"]
  
  tab_Stat_MethCores=rbind.data.frame(tab_Stat_MethCores,tab_Stat_cores)
  
  MethCore_Index=data.frame(average_methylation=-1,MethCores_Index=-1,Groups=NA,Group=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1)
  MethCore_Index$average_methylation=mean(Stat1$Means)
  
  if (dim(Stat_cores)[1]==0) {
    MethCore_Index$MethCores_Index=max(vertex$Freq)
  } else {
    MethCore_Index$MethCores_Index=Stat_cores$Freq
  }
  
  MethCore_Index$Clonality_Index=MethCore_Index$MethCores_Index/MethCore_Index$average_methylation
  MethCore_Index$Group=timing[t]
  MethCore_Index$Groups=unique(Stat1$Tissue[MethCore_Index$Group == Stat1$Group])
  MethCore_Index$Tissue=Groups[g]
  MethCore_Index$Gene=gene
  MethCore_Index$MethCores=unique(my_data3$MethCores)
  quantile=unique(my_data3$quantile)
  quantile=quantile[quantile != ""] 
  
  if (length(quantile !=0)) {
    MethCore_Index$quantile=quantile
  } else { 
    MethCore_Index$quantile=""
  }
  #Replace value grater to 1 with 1
  if (MethCore_Index$Clonality_Index>1) {
    MethCore_Index$Clonality_Index=1
  } else {
    MethCore_Index$Clonality_Index=MethCore_Index$Clonality_Index
  }
  
  MethCore_Index$sd_MethCores_Index=unique(tab_Stat_cores$sd_MethCores_Index)
  MethCore_Index$se_MethCores_Index=unique(tab_Stat_cores$se_MethCores_Index)
  
  MethCore_Index$sd_Clonality_Index=unique(tab_Stat_cores$sd_Clonality_Index)
  MethCore_Index$se_Clonality_Index=unique(tab_Stat_cores$se_Clonality_Index)
  
  MethCores_Index=rbind(MethCores_Index,MethCore_Index)
  
  
  #Legenda Selection
  #legdata1=legend('topright',title = "MethCores Index",legend=c(paste(timing[t],"=",round(MethCore_index$MethCores_index,digits = 2),sep="")),col="black",cex=2)
  #legend(list(x=legdata1$rect$left,y=legdata1$rect$top+0.5),legend=c("Negative Selection, pvalue<=0.01","No Selection","Positive Selection, pvalue<=0.01"),col=c("blue","gray", "red"), lty=1:1, cex=0.8)
  #legend(list(x=legdata1$rect$left,y=legdata1$rect$top+0.5),legend=c("Negative Selection, pvalue<=0.01","No Selection","Positive Selection, pvalue<=0.01"),col=c("blue","gray", "red"), lty=1:1, cex=0.8)
  
  dev.off()

  tab_Stat_MethCores=tab_Stat_MethCores[!is.na(tab_Stat_MethCores$Samples),]
  MethCores_Index=MethCores_Index[!is.na( MethCores_Index$Group),]
  
  #Structure=rbind(Structure,my_data3)
  MethCores_tab1=list(tab_Stat_MethCores=(tab_Stat_MethCores),MethCores_Index=(MethCores_Index),my_data=(my_data1))
  rm(my_data1,my_data2,my_data3,network,vertex,vertex1)
  
  
  return( MethCores_tab1)
}