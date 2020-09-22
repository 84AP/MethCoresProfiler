cores_status_assignment <- function(str1,str2,Complexity,Group,hh,g,tmap1,maxCG_pos,Stat,vertex_cutOFF,dim1,Epi_list1,gene_destination_folder,Groups) {
  
  #Structure=data.frame(From=NA,to=NA,weight=-1,Group=NA,Complexity=NA,Tissue=NA,Evolution=NA)
  #MethCores_Index=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Group=NA,Groups=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1)
  #tab_Stat_MethCores=data.frame(Samples=NA,average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,MethCores=NA,average_MethCores_Index=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,average_Clonality_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1)
  #Description=data.frame(matrix(NA,ncol=(dim(CG_pos)[2]+2)))
  #colnames(Description)=c(CG_pos[1,],"Group","cores")
  
  #Estrai rows with value >0 in tutte le colonne
  #tab_Structure=m1[m1[,c(1:(dim(m1)[2]-1))] >0, ] 
  #tab_Structure= tab_Structure[apply(tab_Structure[,c(1:(dim(m1)[2]-1))],1,function(z) !any(z<0)),] 
  #tab_Structure=tab_Structure[complete.cases(tab_Structure), ]
  #tab_Structure$evolution="Deterministic"
  
  ### Extract cores
  dim3=as.data.frame(str2)
  sel <- apply(dim3[,1:(dim(dim3)[2])],1,function(row) "red" %in% row)
  dim4=dim3[sel,]
  
  # Keep only the lines that are positive in all samples in the group
  #cores=rownames(dim4[rowSums(dim4 == "red")== (length(sample_order)-1), , drop = FALSE])
  dim4$core=rownames(dim4)
  dim5=reshape::melt(dim4,id.vars="core",variable_name="Samples")
  #remove negative rows
  remove=c("blue","white")
  dim5=dim5[!grepl(paste(remove, collapse="|"), dim5$value),]
  
  i <- sapply(dim5, is.factor)
  dim5[i] <- lapply(dim5[i], as.character)
  
  #Check if all samples are inclused
  indx <- unlist(dim5$Samples)
  source=unlist(tmap1$id)
  
  #if  (all(indx %in% source)) {
    #dim5=dim5
  #} else {
   
  # Check if miss some samples
    diff=unique(tmap1$id[!tmap1$id %in% as.character(unique(dim5$Samples))])
    if (length(diff)>=1) {
    for (d in 1:length(diff)) {
    diff1=data.frame(matrix(NA,ncol=3, nrow=1))
    colnames(diff1)=c("core","Samples","value")
    dim1_subset <- as.data.frame(dim1[,diff[d]])
    colnames(dim1_subset)="diff"
    dim1_subset=dim1_subset[order(-dim1_subset$diff), , drop = FALSE]
    
    diff1$Samples=diff[d]
    diff1$core=rownames(dim1_subset)[1]
    diff1$value= "blue"
    #diff1$value= dim1_subset[1,]
      dim5=rbind(dim5,diff1)
    }
    } else {
      dim5=dim5
  }
  # remove Random
    dim5=dim5[!(dim5$Samples=="Random"),]
    
  str_1=as.data.frame(str1)
  str_1$core=rownames(str_1)
  str_1=reshape::melt(str_1,id.vars="core",variable_name="Samples")
  
 dim5 <- merge(dim5,str_1,by=c("core","Samples"),all.x=T)
 colnames(dim5)=c("core","Samples","Selection","Freq")

 if (unique(grepl('_',dim5$Samples))) {
   dim5$Group=sub("_[^_]+$", "",dim5$Samples)
   dim5$Group=gsub('.*_ ?(\\w+)', '\\1', dim5$Group)
 #dim5$Group=(genXtract(dim5$Samples, paste((unique(tmap1$Tissue)[-length(unique(tmap1$Tissue))]),"_",sep=""), "_"))
 } else {
   dim5$Group=dim5$Samples
 }
 dim5$Tissue=sub("_.*","",dim5$Samples)
 #dim5$Tissue <- tmap1$Tissue[match(dim5$Samples, tmap1$id)]
 #### Create a dendrogram of significative epialles
 # Dissimilarity matrix
 df=str1
 # if (max(df)>=vertex_cutOFF) {
 # df=df[ rowSums(df)>= vertex_cutOFF, ]
 # } else {
 #   df=df
 # }
 df1=as.data.frame(t(df))
 df1$Samples=rownames(df1)
 df1$Groups=sub("_[^_]+$", "", df1$Samples)
   
 # before the PCA analysis
 shape=c(15,16,17,3,0,1,2,4)
 shape1=as.data.frame(unique(tmap1$Description))
 colnames(shape1)="Group"
 shape1$shape=shape[1:dim(shape1)[1]]
 shape1$Group=as.factor(shape1$Group)
 shape1$shape=as.integer(shape1$shape)
 #number of duplicates per column
 times=max(unique(sub("M","",tmap1$Rep)))
 shape1=shape1[rep(seq_len(nrow(shape1)),each=times),]
 #shape1=shape1[1:(dim(shape1)[1]-2),]
 shape1$Group=lapply(shape1$Group,as.character)
 shape1$shape=as.numeric(as.character(shape1$shape))

 df1$Samples= factor(df1$Samples, levels = df1$Samples)
 df1$Groups= sub(".*_","",df1$Groups)
 #df1$Groups= as.factor(df1$Groups)
 #df1$Tissue <- tmap1$Tissue[match(rownames(df1), tmap1$id)]
 df1$Tissue=sub("_.*","",df1$Samples)
 df1$shape <- shape1$shape[match(df1$Groups,shape1$Group)]

 res.pca <- prcomp(df1[,1:(dim(df)[1])])
 
 if (Groups[g]=="All") {

palette=rainbow(length(unique(df1$Tissue)))
palette1=as.data.frame(unique(df1$Tissue))
colnames(palette1)="Tissue"
palette1$palette=palette[1:dim(palette1)[1]]

df1$colors <- palette1$palette[match(df1$Tissue,palette1$Tissue)]
 
scores = as.data.frame(res.pca$x) 
scores$Groups=df1$Groups
scores$ColorsTissue=df1$colors
scores$Tissue=df1$Tissue
scores$GroupsShape=df1$shape
scores$GroupsShape[dim(scores)[1]]=11

#ColorsGroups=as.vector(palette1$ColorsGroups)#[as.integer(palette1$Groups)]

var_axis <- apply(res.pca$x, 2, var)

## Proportional variance
var_axis=round(var_axis/sum(var_axis),digits = 2)

summary(res.pca)
PCA1=paste("PC1-Percent variation explained ",(var_axis[1]*100),"%",sep="")
PCA2=paste("PC2-Percent variation explained ",(var_axis[2]*100),"%",sep="")

col_palette <- data.frame(unique(scores)[,c("Tissue","ColorsTissue")])
col_palette =col_palette [order(scores$Tissue),] 

png(paste(gene_destination_folder,"/Groups_Significative_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(scores[,1:2], col=scores$ColorsTissue, pch=scores$GroupsShape,cex=2,
     xlim = c(min(scores$PC1),max(scores$PC1)), ylim=c(min(scores$PC2),max(scores$PC2)), 
     main=paste("PCA ", Groups[g]," ",gene,sep=""),
     xlab=PCA1, ylab=PCA2)

ordiellipse(ord=scores[,1:2], groups = scores$Tissue, kind = "ehull",conf=0.95,draw = "line", col=unique(col_palette$ColorsTissue),#c("blue","red","green","purple"),
           cex=1, label=FALSE,lty = 2,lwd = 1)

#text(scores[,1],scores[,2], labels=rownames(scores), cex= 0.4, pos=3)
# so turn off clipping:

# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")

legend("bottomright",                                # position of legend
       inset=c(-0.315,0),
       legend=unique(scores$Tissue),      # legend display
       pch=21,                                    # point shape
       pt.bg=unique(scores$ColorsTissue),    # point colors
       pt.cex=1.5,                                # point size
       col =unique(scores$ColorsTissue),    # point border color
       title="Tissue")
# 
 legend("topright",                                # position of legend
         inset=c(-0.315,0),
         legend=unique(scores$Groups),      # legend display
         pch=unique(scores$GroupsShape),                                    # point shape
         pt.bg="black",#unique(scores$ColorsTissue),    # point colors
         pt.cex=1.5,                                # point size
         col ="black",#unique(scores$ColorsGroups),    # point border color
         title="Groups")

 #plot(p)
 dev.off()
 
 } else {
   
   palette=rainbow(length(unique(df1$Groups)))
   palette1=as.data.frame(unique(df1$Groups))
   colnames(palette1)="Groups"
   palette1$palette=palette[1:dim(palette1)[1]]
   
   df1$colors <- palette1$palette[match(df1$Groups,palette1$Groups)]
   
   scores = as.data.frame(res.pca$x) 
   scores$Groups=df1$Groups
   scores$ColorsGroups=df1$colors
   scores$Tissue=df1$Tissue
   scores$GroupsShape=df1$shape
   scores$GroupsShape[dim(scores)[1]]=8
   
   #ColorsGroups=as.vector(palette1$ColorsGroups)#[as.integer(palette1$Groups)]
   
   var_axis <- apply(res.pca$x, 2, var)
   
   ## Proportional variance
   var_axis=round(var_axis/sum(var_axis),digits = 2)
   
   summary(res.pca)
   PCA1=paste("PC1-Percent variation explained ",(var_axis[1]*100),"%",sep="")
   PCA2=paste("PC2-Percent variation explained ",(var_axis[2]*100),"%",sep="")
   
   col_palette <- data.frame(unique(scores)[,c("Groups","ColorsGroups")])
   col_palette =col_palette [order(scores$Groups),] 
   
   png(paste(gene_destination_folder,"/Groups_Significative_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
   # Add extra space to right of plot area; change clipping to figure
   par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
   
   plot(scores[,1:2], col=scores$ColorsGroups, pch=scores$GroupsShape,cex=2,
        xlim = c(min(scores$PC1),max(scores$PC1)), ylim=c(min(scores$PC2),max(scores$PC2)), 
        main=paste("PCA ", Groups[g]," ",gene,sep=""),
        xlab=PCA1, ylab=PCA2)
   
   if (length(unique(scores$Groups))==dim(scores)[1]) {
   ordiellipse(ord=scores[,1:2], groups = scores$Groups, kind = "ehull",conf=0.95,draw = "line", col=unique(col_palette$ColorsGroups),#c("blue","red","green","purple"),
               cex=1, label=FALSE,lty = 2,lwd = 1)
   } else {
     ordiellipse(ord=scores[,1:2], groups = scores$Groups, kind = "ehull",conf=0.95,draw = "line", col=unique(col_palette$ColorsGroups),#c("blue","red","green","purple"),
                 cex=1, label=FALSE,lty = 2,lwd = 1)
   }
   
   #text(scores[,1],scores[,2], labels=rownames(scores), cex= 0.4, pos=3)
   # so turn off clipping:
   
   # Add grid lines
   abline(v=0, lty=2, col="grey50")
   abline(h=0, lty=2, col="grey50")
   
   legend("bottomright",                                # position of legend
           inset=c(-0.315,0),
           legend=(unique(scores$Tissue)[-length(unique(scores$Tissue))]),      # legend display
           #pch=c(0,scores[(dim(scores)[1]),17]),                                    # point shape
           #pt.bg=c("white",scores[(dim(scores)[1]),15]),    # point colors
           pt.cex=1.5,                                # point size
           #col =c("white",scores[(dim(scores)[1]),15]),    # point border color
           title="Tissue")
   # 
   legend("topright",                                # sposition of legend
          inset=c(-0.315,0),
          legend=unique(scores$Groups),      # legend display
          pch=unique(scores$GroupsShape),                                    # point shape
          pt.bg=unique(scores$ColorsGroups),    # point colors
          pt.cex=1.5,                                # point size
          col =unique(scores$ColorsGroups),    # point border color
          title="Groups")
   
   #plot(p)
   dev.off()
 }
 
 timing=unique(as.vector(tmap1$Description))[-length(unique(tmap1$Description))]
 
   if (Groups[g]=="All") {
     MethCores_Index2=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Group=NA,Groups=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
     tab_Stat_MethCores2=data.frame(Samples=NA,average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,MethCores=NA,average_MethCores_Index=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,average_Clonality_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
     
     timing1=unique(as.vector(dim5$Tissue)) 
     
        for (t1 in 1:length(timing1)) { 
          for (t in 1:length(timing)) {
       filename=paste(Groups[g],timing1[t1],sep="_")
       dim6=subset(dim5, Group == timing[t] & Tissue == timing1[t1])
       
       if (dim(dim6)[1]==0) next

     MethCores_tab1=MethCores_tab2(tdim5=dim6,timing,t,tmap1,g,filename,Complexity,hh,timing1,t1,Epi_list1,gene_destination_folder,Groups)
     tab_Stat_MethCores1=MethCores_tab1$tab_Stat_MethCores
     MethCores_Index1=MethCores_tab1$MethCores_Index
     
     tab_Stat_MethCores2=rbind(tab_Stat_MethCores2,tab_Stat_MethCores1)
     MethCores_Index2=rbind(MethCores_Index2,MethCores_Index1)
    } 
        } 
     
    # tab_Stat_MethCores=tab_Stat_MethCores2
     #MethCores_Index=MethCores_Index2
     
     } else {
       
       timing1=""
       
       MethCores_Index2=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Group=NA,Groups=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
       tab_Stat_MethCores2=data.frame(Samples=NA,average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,MethCores=NA,average_MethCores_Index=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,average_Clonality_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
       
     for (t in 1:length(timing)) { 
     filename=paste(Groups[g],sep="_")
     MethCores_tab1=MethCores_tab2(tdim5=dim6,timing,t,tmap1,g,filename,Complexity,hh,timing1,t1,Epi_list1,gene_destination_folder)
  
     tab_Stat_MethCores1=MethCores_tab1$tab_Stat_MethCores
     MethCores_Index1=MethCores_tab1$MethCores_Index
     
     tab_Stat_MethCores2=rbind(tab_Stat_MethCores2,tab_Stat_MethCores1)
     MethCores_Index2=rbind(MethCores_Index2,MethCores_Index1)
     
     }
   }
 tab_Stat_MethCores= tab_Stat_MethCores2
 MethCores_Index=MethCores_Index2
 

 tab_Stat_MethCores=tab_Stat_MethCores[!is.na(tab_Stat_MethCores$Samples),]
 write.table(x= tab_Stat_MethCores, file=paste(gene_destination_folder,"/Table_Stat_MethCores_index_",Groups[g],"_",class,"_",gene,".txt", sep=""), row.names=F,col.names=T, quote = F, sep="\t")
 
 MethCores_Index=MethCores_Index[!is.na( MethCores_Index$Group),]
 #MethCores_index$MethCores_index=replace(MethCores_index$MethCores_index, 100 > x >1,1)
 write.table(x= MethCores_Index, file=paste(gene_destination_folder,"/Table_MethCores_index_",Groups[g],"_",class,"_",gene,".txt", sep=""), row.names=F,col.names=T, quote = F, sep="\t")
 
 if (Groups[g]=="All") {
 MethCores_Index$MethCores_Index=round(as.numeric(MethCores_Index$MethCores_Index),digits = 2)
 MethCores_Index$Group <- factor(MethCores_Index$Group, levels = unique(tmap1$Description))
 } else {
   MethCores_Index$Group <- factor(MethCores_Index$Group, levels = MethCores_Index$Group)
 MethCores_Index$MethCores_Index=round(as.numeric(MethCores_Index$MethCores_Index),digits = 2)
}
 # Clonality Index
 png(paste(gene_destination_folder,"/Clonality_Index_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
 
 j1=ggplot(MethCores_Index, aes(x=Group,y=Clonality_Index,fill=Group)) + 
   geom_bar(stat="identity",color="black",position=position_dodge()) +
   geom_errorbar(aes(ymin=Clonality_Index-se_Clonality_Index, ymax=Clonality_Index+se_Clonality_Index), width=.2,position=position_dodge(.9))+
   facet_wrap(~ Groups, scales="free_x") +
   xlab("Samples") +
   ylab("Clonality_Index = % MethCoresFreq / Average Methylation") +
   ylim(0,1.5) + #(as.numeric(1+max(MethCores_Index$se_Clonality_Index))))+
   labs(title =paste("Clonality_Index_",filename,"_",class,"_",gene,sep="")) + 
   labs(fill = "Groups") +
   theme(plot.title = element_text(size = 10))+
   #theme_bw() # Black and white theme
   theme_classic()+ # Classic theme
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
 
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
 
 
 png(paste(gene_destination_folder,"/MethCores_Index_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
  
 j2=ggplot(MethCores_Index, aes(x=Group,y=MethCores_Index,fill=Group)) + 
   geom_bar(stat="identity",color="black",position=position_dodge()) +
   geom_errorbar(aes(ymin=MethCores_Index-se_MethCores_Index, ymax=MethCores_Index+se_MethCores_Index), width=.2,position=position_dodge(.9))+
   facet_wrap(~ Groups, scales="free_x") +
   xlab("Samples") +
   ylab("MethCores_Index = % MethCoresFreq in Total Population") +
   ylim(0,1.5)+
   labs(title =(paste("MethCores_Index_",filename,"_",class,"_",gene,sep=""))) + 
   #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")+
   labs(fill = "Groups") +
   theme(plot.title = element_text(size = 10))+
   #theme_bw() # Black and white theme
   theme_classic()+ # Classic theme
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
 
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

 MethCores=list(MethCores=MethCores_Index,dim5=dim5)
  return(MethCores)
}

