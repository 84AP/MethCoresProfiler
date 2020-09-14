make_graph<- function(tmer,tmer2,sam,summ_stat) {
  #Methylated Grafo
  #Create a df
  tab_mer2=tmer2
  
  #set population name
    class="Methylated"
    class1="Methylation"
    vertex_color="orange"
    n_type="n_Meth"
    CpGs_Freq=ner1
    
  diag(tab_mer2)=0
  yy1 = 0
  #y2= quantile(mer2, 0.7)
  
  #Select the values of the matrix = or> of "y1"
  mcorr1 <-  which(tab_mer2>yy1,arr.ind=TRUE) #extract row and col numbers including those of the diagonal
  zcorr1 <- tab_mer2[tab_mer2>yy1] #extract values associated with the rows and columns of mcorr including those of the diagonal
  cb1=cbind(mcorr1,zcorr1)
  #Create a dataframe with the extracted values and the names of the rows and columns

  tab_df1 <- data.frame(cbind(unlist(cb1), #this creates the desired data frame
                              rownames(tab_mer2)[mcorr1[,1]],
                              colnames(tab_mer2)[mcorr1[,2]]
  ))
  
  #it eliminates the diagonal, the identical positions
  tab_df1<-tab_df1[!(tab_df1$V4==tab_df1$V5),]
  
  #extract the combinations
  tab_list=tab_df1[, 4:5]

  #Merge the columns again
  tab_list3=tab_list
  if (dim(tab_list3)[1]==0) {
    gtab_list3=as.data.frame(t(combn(CG_pos[1,],m = 2)))
    #gtab_list3=gtab_list3[1:10,]
    #Convert V5 column numbers to 3 digits (e.g. 035)
    g1tab_list3=as.data.frame(gtab_list3$V1)
    
    #Convert V5 column numbers to 3 digits (e.g. 035)
    g2tab_list3=as.data.frame(gtab_list3$V2)
    
    tab_list3=gtab_list3
    colnames(tab_list3)[1]="From"
    colnames(tab_list3)[2]="to"
    
    tab_list3$weight=0
    tab_list3$Selection="Selected"
    tab_list3$id=paste(tab_list3$From,tab_list3$to,sep="_")
    crl=tab_list3
  } else {
    
    #Sort the rows in ascending order
    tab_list4=t(apply(tab_list3, 1, sort))
    
    ### Genera Circular graph
    crl=as.data.frame(tab_list4[,c(1:2)])
    crl$weight=tab_df1$zcorr1
    colnames(crl)[1]="From"
    colnames(crl)[2]="to"
    crl$id=paste(crl$From,crl$to,sep="_")
    
    crl=crl[!duplicated(crl[("id")]),]
  }
  #I only take from the 80th percentile onwards
  diag(tab_mer2)=0
  y1 = 0
  #y2= quantile(mer2, 0.7)
  
  #Select the values of the matrix = or> of "y1"
  mcorr <-  which(tab_mer2>y1,arr.ind=TRUE) #extract row and col numbers including those of the diagonal
  zcorr <- tab_mer2[tab_mer2>y1] #extract values associated with the rows and columns of mcorr including those of the diagonal
  cb=cbind(mcorr,zcorr)
  #Create a dataframe with the extracted values and the names of the rows and columns
  
  df <- data.frame(cbind(unlist(cb), #this creates the desired data frame
                         rownames(tab_mer2)[mcorr[,1]],
                         colnames(tab_mer2)[mcorr[,2]]
  ))
  
  #it eliminates the diagonal, the identical positions
  df<-df[!(df$V4==df$V5),]
  
  #extract the combinations
  list=df[, 4:5]
  
  #Merge the columns again
  list3=list
  #Ordina in maniera crescente le righe
  list4=unique(t(apply(list3, 1, sort)))
  list44=(t(apply(list3, 1, sort)))
  
  #write.table(x = list4, file = paste("clust_",sam, ".txt",sep=""), quote=F,sep="\t",row.names=F)
  
  ### Genera Circular graph
  if (dim(list4)[1]==0) {
    glist4=as.data.frame(t(combn(CG_pos[1,],m = 2)))
    
    list4=glist4
    colnames(list4)[1]="From"
    colnames(list4)[2]="to"
    
    list4$weight=0
    list4$Selection="Discarded"
    list4$id=paste(list4$From,list4$to,sep="_")
    crll=list4
    xcrl=merge(crl,crll,by="id",all.x = T)
  } else {
    crll=as.data.frame(list44[,c(1:2)])
    crll$weight=df$zcorr
    colnames(crll)[1]="From"
    colnames(crll)[2]="to"
    crll$id=paste(crll$From,crll$to,sep="_")
    #Sort rows in ascending order
    crll=crll[!duplicated(crll[("id")]),]
    
    i <- sapply( crll, is.factor)
    crll[i] <- lapply(crll[i], as.character)
    crll$weight <- as.numeric(as.character(crll$weight))
    
    ##Highlight selected cores above 70%
    xcrl=merge(crl,crll,by="id",all.x = T)
    i <- sapply( xcrl, is.factor)
    xcrl[i] <- lapply(xcrl[i], as.character)
    
    xcrl$From.y[!is.na(xcrl$From.y)] <- "Selected"
    xcrl$From.y[is.na(xcrl$From.y)] <- "Discarded"
    
    crl=xcrl[,c(2:5)]
    colnames(crl)=c("From","to","weight","Selection")
  }
  #crl$id=NULL
  #Create  all combinations
  combination=as.data.frame(t(combn(CG_pos[1,],m = 2)))
  
  combination1=combination
  
  colnames(combination1)[1]="From"
  colnames(combination1)[2]="to"
  
  combination1$weight=0
  combination1$Selection="Discarded"
  combination1$id=paste(combination1$From,combination1$to,sep="_")

  #Remove existing combinations
  combination2=combination1[!(combination1$id %in% xcrl$id),]
  
  crl=unique(rbind(crl[,1:4],combination2[,1:4]))
  
  #myColors=randomColor(dim(CG_pos)[2],luminosity = "light")
  #### Create a Graph ###
  crl1=data.frame(t(CpGs_Freq[1,c(1:dim(CG_pos)[2])]))
  colnames(crl1)[1]="Freq_CGpos"
  crl1$CGpos=rownames(crl1)
  
  i <- sapply( crl, is.factor)
  crl[i] <- lapply(crl[i], as.character)
  
  crl$weight=as.numeric(crl$weight)
  crl$Freq_From=-1
  crl$Freq_to=-1
  crl$expected_Freq=-1
  crl$observed_Freq=-1
  
  tcrl1=t(crl1)
  colnames(tcrl1)=tcrl1[2,]
  
  #i=1
  #Add Freq to each CG
  for (i in 1:nrow(crl)) {
    fCG=crl[i,]
    CG_From=as.character(fCG$From)
    CG_to=as.character(fCG$to)
    #CG_From1=sub("^0", "", CG_From)
    #select a particular sample from the metamap file
    Freq_From=data.frame(t(subset(tcrl1, select=c(CG_From))))
    Freq_to=data.frame(t(subset(tcrl1, select=c(CG_to))))
    crl$Freq_From[i]=as.numeric(as.character(Freq_From$Freq_CGpos))
    crl$Freq_to[i]=as.numeric(as.character(Freq_to$Freq_CGpos))
  }
  #crl$expected_Freq=((crl$Freq_From)*(crl$Freq_to))/100
  
  #Calcolate observed Freq
  #Create a dataframe to extract the profiles
  cl2=as.data.frame(crl[,1:2])
  cl2[] <- lapply(cl2, as.character)
  #cl2$From=gsub("^0", "", cl2$From)
  #cl2$to=gsub("^0", "", cl2$to)
  
  
  #Search for any combination and apply summary function
  clones=as.data.frame(tmer)
  factor <- sapply(clones, is.factor)
  clones[ factor] <- lapply(clones[ factor], as.character)
  #clones$n_readsCloneForSample <- as.numeric(as.character(clones$n_readsCloneForSample )) 
  
  columns <-c("n_readsCloneForSample","tot_reads_sample","FreqCloneForSample","n","n_readsGroupForSample")
  clones[, columns] <- lapply(columns, function(x) as.numeric(clones[[x]]))
  
  
  #Search for any combination and apply summary function

  summary1= foreach(k=1:nrow(cl2), .combine = rbind) %dopar% {
    #create an empty df
    summary=data.frame(samples=NA, core=NA, Tot=-1, nEpialCore=-1,  FreqEpialCore=-1, nEpialUnici=-1,FreqEpialCoreUnici=-1,FreqEpialMoreFreq=-1)
    #assign name
    summary$samples=sam
    #calculate the number of epialleles
    summary$Tot=sum(clones$n_readsCloneForSample)  
    combs=tidyr::unite_(cl2[k,], paste(colnames(cl2[k,]), collapse="_"), colnames(cl2[k,]),sep="-")
    colnames(combs)="combs"
    #remove duplicates from the list
    combs=vapply(lapply(strsplit(combs$combs, "-"), unique), paste, character(1L), collapse = "-")
    #create a list with the cpg set for each epiallele (you can also put this out of the loop)
    cpg_sets <- strsplit(clones$id_pos,"-")
    #look for indices of sets that include all cpgs of core cl2 [k,] regardless of order
    core_ind <- which(sapply(cpg_sets,function(x) { all(cl2[k,] %in% x)}))
    #extract the core
    core <- clones[core_ind,]
    
    summary=summ_stat(summ = summary, cr = core, cmb = combs) ##### la stringa di ricerca va passata da fuori
  }
  stopImplicitCluster()
  
  crl$expected_Freq=((crl$Freq_From)*(crl$Freq_to))
  #crl$expected_Freq=((crl$Freq_From*crl$Freq_to)/100)
  crl$observed_Freq=(summary1$FreqEpialCore)
  
  crl[is.na(crl)] <- 0
  crl$Tot=b4
  
  #Remove objects
  rm(df)
  
  agaist="expected"
  df2=chi_square(mer1=crl,agaist1=agaist)
  
  crl=cbind(crl,df2[,10:(dim(df2)[2])])
  
  colnames(crl)[10:11] <- c('xSquared', 'pvalue')
  
  #Assign selection type
  for (i in 1:nrow(crl)) {
    ncrl=crl[i,]
    if (ncrl$pvalue<=pvalue & ncrl$observed_Freq>ncrl$expected_Freq) {
      crl$Selection[i]=paste("Positive Selection, pvalue<=",pvalue,sep="")
    } else {
      if (ncrl$pvalue<=pvalue & ncrl$observed_Freq<ncrl$expected_Freq) {
        crl$Selection[i]=paste("Negative Selection, pvalue<=",pvalue,sep="") 
      } else {
        crl$Selection[i]="No Selection"
      }
    }
  }
  
  #Convert character to numeric
  cols.num <- c("weight","Freq_From","Freq_to","expected_Freq","observed_Freq","Tot","xSquared","pvalue")
  crl[cols.num] <- sapply(crl[cols.num],as.numeric)
  
  #Crea oggetto igraph
  g <- graph.data.frame(crl[,1:2], directed=F)
  
  # create a vector of vertex sizes conditional on g*values in df
  # set missing values to size 0 
  V(g)$name=crl1$CGpos
  
  #if (max(crl1$Freq_CGpos)>30){
   # r <- data.frame(g=crl1$CGpos, value=(crl1$Freq_CGpos/3))
  #} else {
    r <- data.frame(g=crl1$CGpos, value=(crl1$Freq_CGpos*10))
    #}
  
  sz <- r$value[match(V(g)$name, r$g)]
  sz[is.na(sz)] <- 0 
  
  #Assign weight to ech vertex
  #if (max(crl$weight)>5){
    E(g)$weight=crl$weight
  #} else {
   # E(g)$weight=crl$weight*10
  #}
  #Color scaling function
  
  #E(g)$color[E(g)$weight==22] <- 'red'
  #Assign colour to ech vertex
  V(g)$name=crl1$CGpos
  #V(g)$color <- "orange"
  la=layout.circle(g)
  
  #create a df  with colors
  color_df=data.frame(matrix(NA, nrow = 3, ncol = 2))
  colnames(color_df)=c("Selection", "Color")
  color_df[,1]=c(paste("Negative Selection, pvalue<=",pvalue,sep=""),"No Selection",paste("Positive Selection, pvalue<=",pvalue,sep=""))
  color_df[,2]=c("blue","gray0","red1")
  
  Selection=unique(crl$Selection)
  #select a particular sample from the metamap file
  mycolors=color_df[color_df$Selection %in% Selection,]
  mycolors1=as.vector(mycolors$Color)
  
  lab.locs <- radian.rescale(x=1:dim(CG_pos)[2], direction=-1, start=0)
  
   #if (max(sz)<=70) {
    #label_dist=2
  #} else {
    label_dist=3
  #}
    if (saveAs=="pdf") {
      pdf(paste("./",destination_folder,"/",sam,"_",class,"_Graph_", gene,".pdf", sep=""),paper = "a4")
    } else { 
  png(paste("./",destination_folder,"/",sam,"_",class,"_Graph_", gene,".png",sep=""),width = 5*150,height = 5*150,res = 200,pointsize = 4)        # smaller font size
    }
      # add the size and border colour vectors to relevant plot arguments     
  plot(g,layout=la,vertex.size=sz,vertex.label=V(g)$name, vertex.color=vertex_color, 
       vertex.frame.color="black", edge.color=as.character(factor(crl$Selection,labels = mycolors1)), edge.width=E(g)$weight,
       vertex.color=V(g)$color, vertex.frame.width=2,rescale=TRUE,vertex.label=TRUE, vertex.label.dist=label_dist,
       vertex.label.degree=lab.locs,vertex.label.cex=1.5, add=FALSE,vertex.label.font=2,margin=0.5)
  title(paste(sam,"_",class,"_Graph_", gene,sep=""),cex.main=2,col.main="black") #
  
  #Legend Methylation
  sizeCut<- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  sizeCutScale <- sizeCut*10
 
  legend('topleft',title = paste(class1," Frequency",sep=""),legend=unique(sizeCut),pt.cex= sizeCutScale,col='black')
  a <- legend('topleft',title = paste(class1," Frequency",sep=""),legend=unique(sizeCut),pt.cex=sizeCutScale/200,col='white', pch=21, pt.bg='white')
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  symbols(x,y,circles=sizeCutScale/200,inches=FALSE,add=TRUE,bg='orange')
  #Legenda methylation
  #legdata=legend('topleft',title = paste(class1," Frequency",sep=""),legend=levels(sizeCut),pt.cex=(scaled/10),col='black',pch=21, pt.bg=vertex_color, cex=1)
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Percent of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  #legdata=legend(list(x=legdata$rect$left,y=legdata$rect$top+0.5),title = "Percent of Methylation", legend=levels(sizeCut),pt.cex=scaled,col='black',pch=21, pt.bg='orange')
  
  #Legenda Selection
  legdata1=legend('bottomright',title = "Types of Selection",legend=c(paste("Negative Selection, pvalue<=",pvalue,sep=""),"No Selection",paste("Positive Selection, pvalue<=",pvalue,sep="")),col=c("blue","gray0","red1"), lty=1, cex=1.2)
  
  #Legenda Widths of Interactions
  #sizeCut<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  sizeCutScale1 <- sizeCut
  
  a=legend('bottomleft',title = paste("Frequency of Interactions",sep=""),legend=unique(sizeCut),lwd=sizeCutScale1,col='black',cex=1)
  #a <- legend('bottomleft',title = paste("Frequency of Interactions",sep=""),legend=unique(sizeCut),lwd=sizeCutScale/200,col='white', cex=1)
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  #symbols(x,y,sizeCutScale/200,inches=FALSE,add=TRUE,bg='orange')
  #legdata2=legend('bottomright', title =paste("Frequency of Interactions",sep=""), legend=levels(sizeCut1),lwd=scaled1, cex=1)
  #legend(list(x=legdata2$rect$left,y=legdata2$rect$top+0.5),title = "Widths of Interactions", legend=levels(sizeCut1),lwd=scaled1)
  #legend(list(x=legdata2$rect$left,y=legdata2$rect$top+0.5), title = "Widths of Interactions", legend=levels(sizeCut1),pt.cex=scaled1,col='black',pch=21, pt.bg='black')
  
  
  dev.off()
  
  crl=data.frame(lapply(crl, as.character), stringsAsFactors=FALSE)
  Positive_crl=crl[grep( "Positive",crl$Selection), ]
  Negative_crl=crl[grep( "Negative",crl$Selection), ]
  
  write.table(x = Positive_crl, file = paste("./",destination_folder,"/",sam,"_Positive_Selection_",class,"_clust_", gene, ".txt",sep=""), quote=F,sep="\t",row.names=F, col.names = T)
  write.table(x = Negative_crl, file = paste("./",destination_folder,"/",sam,"_Negative_Selection_",class,"_clust_", gene, ".txt",sep=""), quote=F,sep="\t",row.names=F, col.names = T)
  write.table(x = crl, file = paste("./",destination_folder,"/",sam,"_All_Links_",class,"_clust_", gene, ".txt",sep=""), quote=F,sep="\t",row.names=F, col.names = T)
}

