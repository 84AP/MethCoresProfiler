epialleles_heatmap2 <- function(Epi_list,Stat,vertex_cutOFF,colorCodes,input_folder2,HighComplexity,Ransom) {
  
  MethCores1=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Group=NA,Groups=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
  MethCores1=MethCores1[0,]
  
  #Extract vector with complexity levels
  
  if (HighComplexity==TRUE) {
  Complexity=complexity(Epi_list)   
  } else {
    Complexity= "dimethyl"
  }
  
  for (hh in 1:length(Complexity))   {
    MethCores=data.frame(average_methylation=-1,MethCores_Index=-1,Clonality_Index=-1,Groups=NA,Group=NA,Tissue=NA,Gene=NA,MethCores=NA,quantile=-1,sd_MethCores_Index=-1,se_MethCores_Index=-1,sd_Clonality_Index=-1,se_Clonality_Index=-1,MoreRepresentativeEpiallele=-1)
    
    Random_summ1=paste("Tab_Epialleli_",Complexity[hh],"_",class,"_",Random_control,"_All_Links_",gene,".txt",sep="")
    
    
    Epi_list1=Epi_list[grep(Complexity[hh], Epi_list, fixed=T)]
    # if Random_simm1 do not exist --> next
    if ((Random_summ1 %in% Random)==FALSE) next
    
    #remove Random_control_Positive
    #Epi_list1=Epi_list1[!str_detect(Epi_list1,pattern="Random_")]
    #Epi_list1=Epi_list1[!str_detect(Epi_list1,pattern="Randomized_M1_All_Links")]
    
    #add Random_control_All_links
    Epi_list1=append(Epi_list1,Random_summ1) 
    
    #Compare Samples cores (Percent)
    em=compare_Samples_epialleles(LIST=Epi_list1,Random_summ1=Random_summ1,map=map,input_folder2=input_folder2)
    #stopCluster(nProcessor)
    
    for (g in 1:length(Groups)) {
      
      if (Groups[g]=="All") {
        tmap1=tmap[1:(dim(tmap)[1]),]
      } else {
        tmap1=tmap[tmap$Group %in% Groups[g],]
        Random1=tmap[tmap$Group %in% "Random$",]
        tmap1=rbind(tmap1,Random1)
      }
      
      #Check if all samples are present in list
      indx <- unlist(tmap1$id)
      source=unlist(colnames(em))
      if (all(indx %in% source)==FALSE) next #all
      #cat(Groups[g])
      
      names.use=c("core","Tot",tmap1$id)
      nm<- em[, names.use]
      #colnames(nm)[dim(nm)[2]]="Random"
      
      #Change column name by match
      #names(nm) <- plyr::mapvalues(names(nm), from = tmap1$SamplesID, to =tmap1$SamplesName)
      
      #Remove rows with zero in all columns
      ##Go through each row and determine if a value is zero
      row_sub = apply(nm, 1, function(row) all(row !=0 ))
      ##Subset as usual
      nm=nm[row_sub,]
      
      m1=nm[,c(3:(dim(nm)[2]))]
      Tot=unique(as.numeric(nm$Tot))
      #m1=m1-m1[,Random_control]
      m1[] <- lapply(m1, function(x) as.numeric(as.character(x)))
      rownames(m1)=nm$core
      
      #Assign colors based on Random control and pvalue
      dim2=Assign_color(m1,Tot,pvalue)
      #dim2=as.matrix(matrix)
      
      write.table(x=dim2, file=paste("./",destination_folder,"/Table_Positive_Epialleles_",Groups[g],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T,col.names=T, quote = F, sep="\t")
      
      #Extract profiles
      cluster=strings_to_BinaryProfiles(df=nm,CG_pos=CG_pos,class=class)
      n= dim(cluster)[2]-dim(CG_pos)[2]
      data=cluster[,c((n+1):dim(cluster)[2])]
      #n2=n1-dim(CG_pos)[2]
      cluster[,c(2:(n))]=sapply(cluster[,c(2:(n))], as.numeric)
      
      data1=data.matrix((cluster[,c(3:(n))]))
      r_names <- cluster$core  # assign labels in column 1 to "rnames"
      rownames(data) <- r_names 
      rownames(data1) <- r_names 
      
      mat_data <- (data[,1:ncol(data)])# transform column 2-5 into a matrix
      
      # Model Based Clustering (Calculate number of clusters)
      #fit <- prof(fit_mclust <- foreach(r = mat_data) %do% try(Mclust(r)))
      fit <- Mclust(mat_data)
      #plot(fit)
      #summary(fit) # display the best model 
      k=as.numeric(fit$G)
      
      ## Plot double heatmap
      dim=data.matrix(mat_data)
      write.table(x=dim, file=paste("./",destination_folder,"/Table_Dendrogram_Epialleles_",Groups[g],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T,col.names=T, quote = F, sep="\t")
      
      #tab_Percento Epiallels
      dim1=data1
      dim1=dim1[match(rownames(dim), rownames(dim1)),]
      write.table(x=dim1, file=paste("./",destination_folder,"/Table_Freq_Epialleles_",Groups[g],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T, col.names=T, quote = F, sep="\t")
      
      kclus <- kmeans(dim, k)
      
      #Create a break and color list
      min=as.numeric(min(dim1))
      max=as.numeric(max(dim1))
      
      #assign intervals
      FunctionIntervalM <- function(min,max) {
        seq(from=0, to =max, length = 5 )
      }
      
      col_breaks=as.numeric(FunctionIntervalM(min,max))
      
      my_palette<- c('white','green','red','cyan','orange','purple','pink','blue','brown','gray','black')
      
      #Generate colors
      col_fun = colorRamp2(col_breaks,colors=my_palette[1:length(col_breaks)], transparency = 0, space = "RGB")
      
      #Associate rows dim1 (Freq epiallels and rige dim2 (cluster red abd blue))
      dim2=dim2[match(rownames(dim1), rownames(dim2)),]
      
      split <- paste0("Cluster\n", kclus$cluster)
      
      #sort colors
      col_dim2=unique(as.vector(unique(dim2)))
      
      #if (length(col_dim2)==3) {
      #  col_dim2=c("blue","red","white")
      #} else if (length(col_dim2)==2) {
      #  col_dim2=c("red","white")
      #}
      col_order=c("blue","red","white")
      
      lgd = Legend(at = c("unmethylated", "methylated"), labels_gp = gpar(fontsize = 5),title = "Clustering of CpG Methylation", type = "points", legend_gp = gpar(col = c("white", "red")),border = TRUE,
                   title_gp = gpar(fontsize = 5, fontface = "bold"))
      lgd1 = Legend(col_fun = col_fun,labels_gp = gpar(fontsize = 5), title = "Frequency of Epiallels", border = TRUE, title_gp = gpar(fontsize = 5, fontface = "bold"))
      lgd2 = Legend(at = c("positive selection", "negative selection","no selection"),title = "Types of Selection",  labels_gp = gpar(fontsize = 5), legend_gp = gpar(fill = c("red","blue","white")),border = TRUE,
                    title_gp = gpar(fontsize = 5, fontface = "bold"))
      pd = packLegend(lgd, lgd1,lgd2,direction = "vertical")
      
      ## Clustering matrix
      row_dend =as.dendrogram(hclust(dist(dim,method="binary")))
      
      if (Groups[g]=="All") {
        #Groups_List=Groups#[-length(Groups)]
        #mdim=list()
        Random_dimm1=data.matrix(dim1[,"Random"])
        colnames(Random_dimm1)="Random"
        
        Random_dimm2=data.matrix(dim2[,"Random"])
        colnames(Random_dimm2)="Random"
        
         dendro1 <- Heatmap(dim,name = "CpG Methylation", show_row_names = FALSE, show_column_names = T, border = TRUE,col = c("white", "red"),
                            column_names_gp = gpar(fontsize =  5, fontface = "bold"),cluster_rows = row_dend,row_dend_width = unit(4, "cm"), cluster_columns = FALSE, show_heatmap_legend = F, row_title = "Cluster")
        
        #for (g2 in 1:length(Groups_List)) {
         # dimm=dim1[,colnames(dim1) %like%  paste(Groups_List[g2],"_",sep="")]
         dimm=dim1
          #dimm=cbind(dimm, Random_dimm1)
          #dimm=dimm[match(rownames(dim), rownames(dimm)),] #paste(Groups_List[g2],sep="")
          dimm1=Heatmap(dimm, cluster_rows = row_dend,name = paste(Groups,sep=""),show_column_names = T,show_row_names = FALSE,border = TRUE,show_heatmap_legend = FALSE,
                        column_names_gp = gpar(fontsize =  5, fontface = "bold"),cluster_columns = FALSE,col = col_fun,column_title=paste(Groups,sep=""),column_title_gp =gpar(fontsize =  5, fontface = "bold"))
          #heatmap_legend_param = list(title = paste("Frequency of Epiallels ",Groups_List[g2],sep=""), title_gp = gpar(fontsize = 5, fontface = "bold"),labels_gp = gpar(fontsize = 5)))
          dendro1=(dendro1+dimm1)
        #}
        
        #for (g3 in 1:length( Groups_List)) {
         # dimmm=dim2[,colnames(dim2) %like%  paste(Groups_List[g3],"_",sep="")]
          dimmm=dim2
          #dimmm=cbind(dimmm, Random_dimm2)
          #dimmm=dimmm[match(rownames(dim), rownames(dimmm)),]
          dimm2=Heatmap(dimmm,cluster_rows = row_dend, name = paste(Groups,sep=""),show_row_names = FALSE, show_column_names = T,border = TRUE,show_heatmap_legend = FALSE,
                        column_names_gp = gpar(fontsize =  5, fontface = "bold"), cluster_columns = FALSE,col=unique(c(dimmm)),column_title=paste(Groups,sep=""), column_title_gp =gpar(fontsize =  5, fontface = "bold")) #col = col_fun1,
          #heatmap_legend_param = list(title = paste("Types of Selection ",Groups_List[g3],sep=""),title_gp = gpar(fontsize = 5, fontface = "bold"),
          #                           labels_gp = gpar(fontsize = 5)))
          dendro1=(dendro1+dimm2)
        #}
        
        ##### Save separately
        dendr_o <- Heatmap(dim, cluster_rows = row_dend, name = "CpG Methylation",width = unit(1, "cm"), show_row_names = FALSE, show_column_names = T, border = TRUE,col = c("white", "red"),
                           column_names_gp = gpar(fontsize =  5, fontface = "bold"),row_dend_width = unit(1, "cm"),cluster_columns = FALSE, show_heatmap_legend = F, row_title = "Cluster")
        
        png(paste("./",destination_folder,"/Clustering_CpGs_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
            width = 5*250,        # 5 x 300 pixels
            height = 5*250,
            res = 300,            # 300 pixels per inch
            pointsize = 4)        # smaller font size
        
        pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
        #grid.text(paste("Clusters_Epialleles_",Group[g],"_",class,"_",Complexity[hh],"_",gene1[1],sep=""), vp = viewport(height = 0, layout.pos.col = 1))
        draw(dendr_o, newpage = FALSE, annotation_legend_list = lgd)
        #column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
        upViewport()
        dev.off()
        
        dendr_o1= NULL
        #for (g_2 in 1:length( Groups_List)) {
          #dim_m=dim1[,colnames(dim1) %like%  paste(Groups_List[g_2],"_",sep="")]
        dim_m=dim1  
        #dim_m=cbind(dim_m, Random_dimm1)
          #dim_m=dim_m[match(rownames(dim), rownames(dim_m)),]
          dim_m1=Heatmap(dim_m,cluster_rows = row_dend, name = Groups,show_column_names = T,show_row_names = FALSE,border = TRUE,show_heatmap_legend = FALSE,show_row_dend = FALSE,
                         column_names_gp = gpar(fontsize =  5, fontface = "bold"),cluster_columns = FALSE,col = col_fun,column_title=Groups,column_title_gp =gpar(fontsize =  5, fontface = "bold"))
          #heatmap_legend_param = list(title = paste("Frequency of Epiallels ",Groups_List[g2],sep=""), title_gp = gpar(fontsize = 5, fontface = "bold"),labels_gp = gpar(fontsize = 5)))
          dendr_o1=(dendr_o1+dim_m1)
        #}
        png(paste("./",destination_folder,"/Clustering_EpiallelesFreq_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
            width = 5*250,        # 5 x 300 pixels
            height = 5*250,
            res = 300,            # 300 pixels per inch
            pointsize = 4)        # smaller font size
        
        pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
        #grid.text(paste("Clusters_Epialleles_",Group[g],"_",class,"_",Complexity[hh],"_",gene1[1],sep=""), vp = viewport(height = 0, layout.pos.col = 1))
        draw(dendr_o1, newpage = FALSE, annotation_legend_list = lgd1)
        #column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
        upViewport()
        dev.off()
        
        dendr_o2= NULL
        #for (g_3 in 1:length( Groups_List)) {
         # dimm_m=dim2[,colnames(dim2) %like%  paste(Groups_List[g_3],"_",sep="")]
        dimm_m=dim2
          #dimm_m=cbind(dimm_m, Random_dimm2)
          col1m=unique(as.vector(unique(dimm_m)))
          #col1[sort(order(y)[x])]
          col1m=col1m[order(match(col1m,col_order))]
          #Associa righe dim1 (Freq epiallels e rige dim2(cluster red abd blue))
          #dimm_m=dimm_m[match(rownames(dim), rownames(dimm_m)),] #col=col_dim2
          dim_m2=Heatmap(dimm_m,cluster_rows = row_dend, name = Groups,show_row_names = FALSE, show_column_names = T,border = TRUE,show_heatmap_legend = FALSE,show_row_dend = FALSE,
                         column_names_gp = gpar(fontsize =  5, fontface = "bold"), cluster_columns = FALSE,col=col1m,column_title=Groups, column_title_gp =gpar(fontsize =  5, fontface = "bold")) #col = col_fun1,
          #heatmap_legend_param = list(title = paste("Types of Selection ",Groups_List[g3],sep=""),title_gp = gpar(fontsize = 5, fontface = "bold"),
          #                           labels_gp = gpar(fontsize = 5)))
          dendr_o2=(dendr_o2+dim_m2)
        #}
        
        png(paste("./",destination_folder,"/Clustering_Selection_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
            width = 5*250,        # 5 x 300 pixels
            height = 5*250,
            res = 300,            # 300 pixels per inch
            pointsize = 4)        # smaller font size
        
        pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
        #grid.text(paste("Clusters_Epialleles_",Group[g],"_",class,"_",Complexity[hh],"_",gene1[1],sep=""), vp = viewport(height = 0, layout.pos.col = 1))
        draw(dendr_o2, newpage = FALSE, annotation_legend_list = lgd2)
        #column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
        upViewport()
        dev.off()
        
      } else {
        col1=unique(as.vector(unique(dim2)))
        #col1[sort(order(y)[x])]
        col1=col1[order(match(col1,col_order))] #col=col_dim2
        dendro1 <- Heatmap(dim, cluster_rows = row_dend, name = "CpG Methylation",width = unit(1, "cm"), show_row_names = FALSE, show_column_names = T, border = TRUE,col = c("white", "red"),
                           column_names_gp = gpar(fontsize =  5, fontface = "bold"),row_dend_width = unit(1, "cm"),cluster_columns = FALSE, show_heatmap_legend = F, row_title = "Cluster")+
          Heatmap(dim1,cluster_rows = row_dend, name = "Frequency of Epiallels",show_column_names = T,show_row_names = FALSE,border = TRUE,show_heatmap_legend = FALSE,show_row_dend = FALSE,
                  column_names_gp = gpar(fontsize =  5, fontface = "bold"),cluster_columns = FALSE,col = col_fun,column_title="Epialleles Freq",column_title_gp =gpar(fontsize =  5, fontface = "bold"))+
          #heatmap_legend_param = list(title = "Frequency of Epiallels", title_gp = gpar(fontsize = 5, fontface = "bold"),labels_gp = gpar(fontsize = 5)))+
          Heatmap(dim2,cluster_rows = row_dend, name = "Types of Selection",show_row_names = FALSE, show_column_names = T,border = TRUE,show_heatmap_legend = FALSE,show_row_dend = FALSE,
                  column_names_gp = gpar(fontsize =  5, fontface = "bold"), cluster_columns = FALSE,col=col1,column_title="chi_square Selection", column_title_gp =gpar(fontsize =  5, fontface = "bold")) #col = col_fun1,
        #  heatmap_legend_param = list(title = "Types of Selection",title_gp = gpar(fontsize = 5, fontface = "bold"),
        #                             labels_gp = gpar(fontsize = 5)))
      
      png(paste("./",destination_folder,"/Clusters_Epialleles_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
          width = 5*250,        # 5 x 300 pixels
          height = 5*250,
          res = 300,            # 300 pixels per inch
          pointsize = 4)        # smaller font size
      
      pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
      #grid.text(paste("Clusters_Epialleles_",Group[g],"_",class,"_",Complexity[hh],"_",gene1[1],sep=""), vp = viewport(height = 0, layout.pos.col = 1))
      draw(dendro1, newpage = FALSE, annotation_legend_list = list(pd))
      #column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
      upViewport()
      dev.off()
    }

    
      if (Complexity[hh]=="dimethyl") {
        MethCores=cores_status_assignment(str1=dim1,str2=dim2,Complexity=Complexity,hh=hh,g=g,Group=Group,tmap1=tmap1,maxCG_pos=maxCG_pos,
                                          Stat=Stat,vertex_cutOFF=vertex_cutOFF,dim1=dim1,Epi_list1=Epi_list1) 
        
        
    ##Dendrogram only for if Group==ALL
        if (Groups[g]=="All") {
          ### Dendrogram with Phylogenetic
          dendrogram =as.data.frame(MethCores$MethCores)#MethCores$dim5 or MethCores$MethCores
          #check if timing = for each group
          #a <- dendrogram %>% group_by(Group) %>% mutate(n_tissues = length(unique(Groups))) %>% ungroup()
          dendrogram = dendrogram %>% 
            dplyr::group_by(Group) %>% 
            dplyr::mutate(n_tissues = n_distinct(Groups)) %>% 
            ungroup()
          
          ##Remove Timing diffenrenti
          dendrogram1<-dendrogram[!(dendrogram$n_tissues<length(unique(dendrogram$Groups))),]
          
          dendrogram1$Samples=paste(dendrogram1$Groups,dendrogram1$Group,sep="_")
          
          # check how many "_" there are in sample name
          for (nnn in 1:nrow(dendrogram1)) {
          if(str_count(dendrogram1$Samples[nnn],"\\_")>1) {
            SAM=strsplit(dendrogram1$Samples[nnn], "\\_")[[1]]
            dendrogram1$Samples[nnn]=paste(SAM[1],SAM[length(SAM)],sep="_")
          } else {
            dendrogram1$Samples[nnn]=dendrogram1$Samples[nnn]
          }
          }
          #dendrogram1=dendrogram[,c(1:2,4)]
          ## If MethCores$MethCores
          dendrogram1=dendrogram1[,c(4,5,8,16)]
          colnames(dendrogram1)=c("Timing","Groups","core","Samples")
          #dendrogram1=dendrogram[,c(1:2,4)]# +Freq
          
          #dendrogram1[,"Freq"]= dendrogram1[,"Freq"]*100
          
          #dendrogram1[,"Freq"]= dendrogram1[,"Freq"]*100
          #dendrogram2 <- cast(dendrogram1, core~Samples, mean,value = 'Freq')
          #dendrogram2[is.na(dendrogram2)] <- 0
          #dendrogram2[,"Random"]=0
          dendrogram2=strings_to_BinaryProfiles(df=dendrogram1,CG_pos=CG_pos,class=class)
          dendrogram2$core=NULL
          dendrogram3=dendrogram2
          #dendrogram3=dendrogram2[,-c(1:2)]
          rownames(dendrogram3)=dendrogram3$Samples
          
          ### Calcola numero di cluster
          #fit1 <- Mclust(dendrogram3)
          #plot(fit)
          #summary(fit) # display the best model 
          #k1=as.numeric(fit1$G)
          
          #kclus1 <- kmeans(dim1, k1)
          
          dendrogram3$strings=apply( dendrogram3[ ,4:(dim(dendrogram3)[2]) ] , 1 , paste , collapse = "" )
          rownames(dendrogram3)=dendrogram3$Samples
          labels=rownames(dendrogram3)
          strings =dendrogram3$strings
          
          # compute hamming distance
          dm <- stringdistmatrix(strings, strings, method = "hamming")
          # convert to distance object, add labels
          colnames(dm) <- labels
          rownames(dm) <- labels
          #colorCodes <- c(CB="red", CX="green", HIPP="yellow3")
          
          d <- dist(dm)
          # hierarchical clustering
          hcl <- hclust(d)
          # plot the tree
          #phylo <- as.phylo(hcl,)
          #phylo$tip.label <- hcl$labels
          if (length(d)>1) {
           png(paste("./",destination_folder,"/Phylogenetic_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
              width = 5*250,        # 5 x 300 pixels
              height = 5*250,
              res = 500,            # 300 pixels per inch
              pointsize = 4,units = c("px"))        # smaller font size
          par(mar=c(9,9,9,9))
        #plot(phylo,main = paste("Phylogenetic_",Groups[g],"_",select,"_",class,"_",Complexity[hh],"_",gene,".txt",sep=""))
        plot(as.phylo(hcl), main = paste("./",destination_folder,"/Phylogenetic_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".txt",sep=""),tip.color=colorCodes[sub("_.*","",rownames(dm))], type="cladogram",
             cex = 1 ,label.offset = 0.2) 
        axisPhylo(1, las = 1,backward=TRUE)
        dev.off()
          }
         ##Label order
         prefix=hcl$labels[sort(order(hcl$labels)[hcl$order])]

        ##Write phylogram
         #k3=cutree(phylo, h = max(hcl$height))
         write.hclust(hcl, file =  paste("./",destination_folder,"/Phylogenetic_Cluster_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".txt",sep=""), prefix = prefix,h=max(hcl$height)) 
         write.table(x=dm, file=paste("./",destination_folder,"/Table_Phylogenetic_Cluster_",Groups[g],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T,col.names=T, quote = F, sep="\t")
         
          
          ### Dendrogram with all Freq
         nm1=as.data.frame(t(nm))
         colnames(nm1)=as.character(nm1[1,])
         nm1=nm1[-c(1:2),]
         
         # check how many "_" there are in sample name
         for (nnnn in 1:nrow(nm1)) {
           if(str_count(rownames(nm1)[nnnn],"\\_")>2) {
             SAM1=strsplit(rownames(nm1)[nnnn], "\\_")[[1]]
             rownames(nm1)[nnnn]=paste(SAM1[1],SAM1[3],SAM1[4],sep="_")
           } else {
             rownames(nm1)[nnnn]=rownames(nm1)[nnnn]
           }
         }
         
         nm1$Samples=rownames(nm1)
         nm1$Timing=sub(".*_P","",nm1$Samples)
         nm1$Timing=paste("p",sub("_.*","",nm1$Timing),sep="")
         nm1$Groups=sub("_.*","",nm1$Samples)
        
         # Compute Euclidean distance between samples
         # Perfor clustering with hclust
         hc=hclust(dist(nm1[ ,(dim(nm1)[2]-3)] , diag=TRUE))
         hc$labels=rownames(nm1)
         dhc <- as.dendrogram(hc)

         #groupCodes <- c(rep("CB",12), rep("CX",12), rep("HIPP",12))
         #rownames(sample) <- make.unique(groupCodes)
         #colorCodes <- c(CB="red", CX="green", HIPP="yellow3")
         # Assigning the labels of dendrogram object with new colors:
         #labels_colors(dhc) <- colorCodes[groupCodes][order.dendrogram(dhc)]
         labels_colors(dhc) <- colorCodes[sub("_.*","",rownames(nm1))][order.dendrogram(dhc)]

         if (length(hc)>1) {
         # And the plot
         png(paste("./",destination_folder,"/Dendrogram_",Groups[g],"_FreqEpiallels_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
             width = 5*250,        # 5 x 300 pixels
             height = 5*250,
             res = 500,            # 300 pixels per inch
             pointsize = 4,units = c("px"))        # smaller font size
         par(mar=c(9,9,9,9))
         plot(dhc,horiz = TRUE ,main = paste("Dendrogram_",Groups[g],"_FreqEpiallels_",class,"_",Complexity[hh],"_",gene,sep=""),
              cex = 1)
       dev.off()
         }
         
         ### Dendrogram with Freq of significant epiallles 
        #core=MethCores$MethCores$MethCores[1]
         dendrogram =as.data.frame(MethCores$MethCores)#MethCores$dim5 or MethCores$MethCores
         #check if timing = for each group
         #a <- dendrogram %>% group_by(Group) %>% mutate(n_tissues = length(unique(Groups))) %>% ungroup()
         dendrogram = dendrogram %>% 
           dplyr::group_by(Group) %>% 
           dplyr::mutate(n_tissues = n_distinct(Groups)) %>% 
           ungroup()
         
         ##Remove Timing diffenrenti
         dendrogram1<-dendrogram[!(dendrogram$n_tissues<length(unique(dendrogram$Groups))),]
         
         dendrogram1$Samples=paste(dendrogram1$Groups,dendrogram1$Group,sep="_")
         # check how many "_" there are in sample name
         for (nnnnn in 1:nrow(dendrogram1)) {
           if(str_count(dendrogram1$Samples[nnnnn],"\\_")>1) {
             SAM2=strsplit(dendrogram1$Samples[nnnnn], "\\_")[[1]]
             dendrogram1$Samples[nnnnn]=paste(SAM2[1],SAM2[length(SAM2)],sep="_")
           } else {
             dendrogram1$Samples[nnnnn]=dendrogram1$Samples[nnnnn]
           }
         }
         #dendrogram1=dendrogram[,c(1:2,4)]
        #dendrogram1=dendrogram[,c(1:2,4)]
        ## If MethCores$MethCores
        dendrogram1=dendrogram1[,c(4,5,8,2,16)]
        colnames(dendrogram1)=c("Timing","Groups","core","Freq","Samples")
        #dendrogram1=dendrogram[,c(1:2,4)]# +Freq

        dendrogram1[,"Freq"]= dendrogram1[,"Freq"]*100
        dendrogram2 <- cast(dendrogram1, core~Samples, mean,value = 'Freq')
        dendrogram2[is.na(dendrogram2)] <- 0
        #dendrogram2[,"Random"]=0
        #dendrogram2=strings_to_BinaryProfiles(df=dendrogram1,CG_pos=CG_pos,class=class)
 
        #dendrogram3=dendrogram2[rep(row.names(dendrogram2), dendrogram2$Freq), 3:dim(dendrogram2)[2]]
        #dendrogram3=t(dendrogram2[,c(2:(dim(dendrogram2)[2]))])
        #dendrogram3=(dendrogram2[,c(3:(dim(dendrogram2)[2]))])
        #dendrogram3=dendrogram2
        dendrogram3=dendrogram2[,-1]
        #rownames(dendrogram3)=paste(dendrogram2$Samples,1:dim(dendrogram2)[1],sep="__")
        dendrogram4=as.data.frame(dendrogram3)
        mat=t(dendrogram4)
        #mat = mtcars
        d = (1 - vegdist(mat, method="bray")) * 100
        #colorCodes <- c(CB="red", CX="green", HIPP="yellow3")
        h = hclust(d)
        hc <- as.dendrogram(h)
        
        labels_colors(hc) <- colorCodes[sub("_.*","",rownames(mat))][order.dendrogram(hc)]

        if (length(d)>1) {
        png(paste("./",destination_folder,"/Dendrogram_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".png",sep=""),
            width = 5*250,        # 5 x 300 pixels
            height = 5*250,
            res = 500,            # 300 pixels per inch
            pointsize = 4,units = c("px"))        # smaller font size
          par(mar=c(9,9,9,9))
        plot(hc, horiz = TRUE ,main = paste("Dendrogram_",Groups[g],"_SignificantEpialleles_",select,"_",class,"_",Complexity[hh],"_",gene,sep=""),#, hang = -1)
             cex = 1)
             dev.off()
        
        ##Write dendrogram
        ##Label order
        prefix=h$labels[sort(order(h$labels)[h$order])]
        write.hclust(h, file =  paste("./",destination_folder,"/Dendrogram_Cluster_",Groups[g],"_",class,"_",Complexity[hh],"_",gene,".txt",sep=""), prefix = prefix,h=max(hcl$height)) 
        write.table(x=mat, file=paste("./",destination_folder,"/Table_Dendrogram_Cluster_",Groups[g],"_",class,"_",Complexity[hh],"_",gene, ".txt", sep=""), row.names=T,col.names=T, quote = F, sep="\t")
        }
        }
        
        
        #####Extract CpGs from MethCores
        MethCores=MethCores$MethCores
        #CpGs=unique(unlist(strsplit(MethCores$MethCores, "-")))
        CpGs=as.character(unlist(CG_pos[1,]))
        #Create a empty df
        tab_MethCores=data.frame(CpGs=NA,Groups=NA,MethCores_Index=-1,id=NA)
        tab_MethCores=tab_MethCores[0,]
        
        #MethCores=MethCores$MethCores
        for (c in 1:nrow(MethCores)) {
          tab=MethCores[c,]
          tab_cores=unlist(strsplit(tab$MethCores, "-"))
          
          tab_Meth=data.frame(matrix(NA,ncol=3,nrow=length(tab_cores)))
          colnames(tab_Meth)=c("CpGs","Groups","MethCores_Index")
          
          #Create a empty df
          # MethCores_Index=percent of cores in total population
          tab_Meth$CpGs=tab_cores
          tab_Meth$MethCores_Index=tab$MethCores_Index
          #tab_Meth$id=paste(tab$Groups,tab$Group,sep="_")
          
          
          if (length(tab_cores)==length(CpGs)) {
            tab_Meth=tab_Meth 
            tab_Meth$Groups=tab$Group
            if (unique(grepl('_',tmap1$SamplesID)[-length(unique(tmap1$SamplesID))])) {
              tab_Meth$id=paste(tab$Groups,tab$Group,sep="_")
            } else {
              tab_Meth$id=tab$Group
            }
            tab_MethCores=rbind(tab_MethCores,tab_Meth)
          } else {
            diff_CpGs=setdiff(CpGs,tab_cores)
            tab_Meth[nrow(tab_Meth) + (length(diff_CpGs)),] = NA
            tab_Meth$CpGs[(length(tab_cores)+1):(dim(tab_Meth)[1])]=diff_CpGs
            tab_Meth$MethCores_Index[(length(tab_cores)+1):(dim(tab_Meth)[1])]=0
            tab_Meth$Groups=tab$Group
            if (unique(grepl('_',tmap1$SamplesID)[-length(unique(tmap1$SamplesID))])) {
              tab_Meth$id=paste(tab$Groups,tab$Group,sep="_")
            } else {
              tab_Meth$id=tab$Group
            }
            
            tab_MethCores=rbind(tab_MethCores,tab_Meth)
            
          }
          
        }
        
        cores=CpGs
        Stat$Group=as.character(Stat$Tissue)
        
        if (Groups[g]=="All") {
          gStat=Stat
          gStat=gStat[!grepl("Random$", gStat$Groups),]
        } else {
          gStat=Stat[Stat$Tissue %in% Groups[g],]
        }
        gStat$Family=gStat$Tissue
        gStat$Group=sub("synthetic_","",gStat$Samples)
        gStat$Group=sub("_.*","",gStat$Group)
        #gStat$id=sub("_[^_]+$","",gStat$Samples)
        aggregateGroups=aggregate(gStat[,1:(dim(CG_pos)[2])], list(gStat$Groups), mean)
        #aggregateGroups=aggregateGroups[!grepl("Random", aggregateGroups$Group.1),]
        #aggregateGroups$Group.1=sub(".*?_","",aggregateGroups$Group.1)
        aggregateGroups_cores <- aggregateGroups[,c("Group.1", cores)]
        aggregateGroups_cores=reshape::melt(aggregateGroups_cores, id="Group.1")
        colnames(aggregateGroups_cores)=c("Groups","CpGs","Freq")
        aggregateGroups_cores$id=aggregateGroups_cores$Groups
        aggregateGroups_cores$Family=sub("_.*","", aggregateGroups_cores$Groups)
        aggregateGroups_cores$Groups=sub(".*_","", aggregateGroups_cores$Groups)
        aggregateGroups_cores$id=paste(aggregateGroups_cores$Family,aggregateGroups_cores$Groups,sep="_")
        #aggregateGroups_cores$id=paste( aggregateGroups_cores$Family, aggregateGroups_cores$Groups,sep="_")
        tab_MethCores=merge(tab_MethCores,aggregateGroups_cores,by=c("Groups","CpGs","id"),all=T)
        
        # Variation of CpG of cores respect to percent of CpG methylation
        tab_MethCores$CpG_Variation=tab_MethCores$MethCores_Index/tab_MethCores$Freq
        #tab_MethCores$Family=sub("_.*","",tab_MethCores$Groups)
        
        tab_MethCores[is.na(tab_MethCores)]=0
        
        #Replace value grater to 1 with 1
        tab_MethCores$CpG_Variation[tab_MethCores$CpG_Variation>1] <-1
        
        
        
        write.table(x=  tab_MethCores, file=paste("./",destination_folder,"/Table_CpG_Contribution_",Groups[g],"_",class,"_",gene,".txt", sep=""), row.names=F,col.names=T, quote = F, sep="\t")
        
        #for (f in 1:length(Family)) {
        
        #tab_MethCores1=tab_MethCores[tab_MethCores$Family %like% Family[f]]
        yearly_counts <- tab_MethCores %>% group_by(CpGs, Groups,CpG_Variation,Family,id) %>% tally
        colnames(yearly_counts)[2]="GGroups"
        
        #if (max(yearly_counts$n)==1) next
        
        #if (length(unique(yearly_counts$Groups))==length(unique(yearly_counts$Family))) {
        if (overtime==TRUE) { 
          X=factor(yearly_counts$Groups,levels = unique(yearly_counts$Groups))
          X1="GGroups"
        } else {
          X=factor(yearly_counts$CpGs,levels = unique(yearly_counts$CpGs))
          X1="CpGs"
          #X=factor(yearly_counts$Groups,levels = unique(tmap1$Description))
          #X1="Groups"
        }
        
        Family=as.character(unique( tab_MethCores$Family))
        GGroups=as.character(unique( yearly_counts$GGroups))
        
        if ((type_Tissue=="same")==TRUE) { 
          #GGroups=unique(yearly_counts$Groups)
          
          png(paste("./",destination_folder,"/MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 250,pointsize = 8) 
          
          j2=ggplot(yearly_counts,aes(X,CpG_Variation, group=CpGs)) + geom_point(aes(colour=CpGs)) + geom_line(aes(colour=CpGs))+
            facet_wrap(~ GGroups, scales="free_x") + #~ Family
            ylim(0,1)+
            xlab(X1) +
            ylab("Freq MethCores_CpG / Freq Methylated CpGs") +
            ggtitle(paste("MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,sep="")) + 
            scale_color_discrete(breaks=cores)+
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
          
          #if ((type_Tissue=="same")==TRUE) { 
          # #GGroups=unique(yearly_counts$Groups)
          #  j2 +facet_wrap(~ Family, scales="free_x") #~ Family
          #} else {
          # j2 +facet_wrap(~ GGroups, scales="free_x") #~ Family
          #}
          
          #j2 +facet_wrap(~ Family, scales="free_x") + #~ Family
          
          
          addSmallLegend <- function(myPlot, pointSize = 6, textSize = 5, spaceLegend = 0.1) {
            myPlot +
              guides(shape = guide_legend(override.aes = list(size = pointSize)),
                     color = guide_legend(override.aes = list(size = pointSize))) +
              theme(legend.title = element_text(size = textSize), 
                    legend.text  = element_text(size = textSize),
                    legend.key.size = unit(spaceLegend, "lines"))
          }
          j2=addSmallLegend(j2)
          plot(j2)
          dev.off()
          
        } else {
          
          png(paste("./",destination_folder,"/MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 250,pointsize = 8) 
          
          j2=ggplot(yearly_counts,aes(X,CpG_Variation, group=CpGs)) + geom_point(aes(colour=CpGs)) + geom_line(aes(colour=CpGs))+
            facet_wrap(~ Family, scales="free_x") + #~ Family
            ylim(0,1)+
            xlab(X1) +
            ylab("Freq MethCores_CpG / Freq Methylated CpGs") +
            ggtitle(paste("MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,sep="")) + 
            scale_color_discrete(breaks=cores)+
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
          
          #if ((type_Tissue=="same")==TRUE) { 
          # #GGroups=unique(yearly_counts$Groups)
          #  j2 +facet_wrap(~ Family, scales="free_x") #~ Family
          #} else {
          # j2 +facet_wrap(~ GGroups, scales="free_x") #~ Family
          #}
          
          #j2 +facet_wrap(~ Family, scales="free_x") + #~ Family
          
          
          addSmallLegend <- function(myPlot, pointSize = 6, textSize = 5, spaceLegend = 0.1) {
            myPlot +
              guides(shape = guide_legend(override.aes = list(size = pointSize)),
                     color = guide_legend(override.aes = list(size = pointSize))) +
              theme(legend.title = element_text(size = textSize), 
                    legend.text  = element_text(size = textSize),
                    legend.key.size = unit(spaceLegend, "lines"))
          }
          j2=addSmallLegend(j2)
          plot(j2)
          dev.off()
          
        }
        
        # png(paste("./",destination_folder,"/MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,".png",sep=""), width = 5*300,height = 5*300,res = 250,pointsize = 8) 
        # 
        # j2=ggplot(yearly_counts,aes(X,CpG_Variation, group=CpGs)) + geom_point(aes(colour=CpGs)) + geom_line(aes(colour=CpGs))+
        #   facet_wrap(~ Family, scales="free_x") + #~ Family
        #   ylim(0,1)+
        #   xlab(X1) +
        #   ylab("Freq MethCores_CpG / Freq Methylated CpGs") +
        #   ggtitle(paste("MethCores_Entanglement_plot_",Groups[g],"_",class,"_",gene,sep="")) + 
        #   scale_color_discrete(breaks=cores)+
        #   theme_classic()+
        #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6,color="darkred"), legend.position = "bottom", legend.box = "horizontal")
        # 
        # #if ((type_Tissue=="same")==TRUE) { 
        #  # #GGroups=unique(yearly_counts$Groups)
        # #  j2 +facet_wrap(~ Family, scales="free_x") #~ Family
        # #} else {
        #  # j2 +facet_wrap(~ GGroups, scales="free_x") #~ Family
        # #}
        # 
        # #j2 +facet_wrap(~ Family, scales="free_x") + #~ Family
        #   
        #   
        # addSmallLegend <- function(myPlot, pointSize = 6, textSize = 5, spaceLegend = 0.1) {
        #   myPlot +
        #     guides(shape = guide_legend(override.aes = list(size = pointSize)),
        #            color = guide_legend(override.aes = list(size = pointSize))) +
        #     theme(legend.title = element_text(size = textSize), 
        #           legend.text  = element_text(size = textSize),
        #           legend.key.size = unit(spaceLegend, "lines"))
        # }
        # j2=addSmallLegend(j2)
        # plot(j2)
        # dev.off()
        # 
        
        rm(list = ls()[!ls() %in% c("HighComplexity","input_folder1","input_folder2","type_Tissue","destination_folder","FOLDERS","colorCodes","destination_folder","overtime","Epi_list1","population_weight","MethCores_tab2","cores_cutOff","MethCoresIndex","MethCores_tab","Assign_color","MethCores1","MethCores",
                                    "maxCG_pos","hh","em","Epi_list","CG_pos","Complex","Complex1","Groups","tmap","Random_control","Selection","class","Complexity","pvalue",
                                    "sample_order", "map", "gene1", "strings_to_BinaryProfiles","list","list1","Mlist","Mplist","gene","Random_summ1","Random_control",
                                    "plist","plist1","compare_Samples_cores","compare_Samples_epialleles","Complexity","complexity","Stat","vertex_cutOFF","nProcessor",
                                    "cores_heatmap","epialleles_heatmap","cores_status_assignment","minimal_common_structure2","Structure1","LIST","MethCores1","MethCores")])    
        
      } else {
        MethCores=MethCores
      }
      
      MethCores1=rbind(MethCores1,MethCores)
      
    }
    
    rm(em,MethCores)
  }
  
  return(MethCores1)
}

