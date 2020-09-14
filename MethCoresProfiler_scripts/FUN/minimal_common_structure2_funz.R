minimal_common_structure2 <- function(cor,tmap1,timing,t,v) {
  
  # tab_cor=data.frame(core=NA,Samples=NA,Evolution=NA)
  # ccor=unique(cor$core)
  # 
  # for ( c in 1:length(ccor)) {
  #   t_cor=data.frame(core=NA,Samples=NA,Evolution=NA)
  #   cor1=subset(cor, core ==ccor[c])
  #   cor1$Samples <- as.character(cor1$Samples)
  #   corSamples=subset(tmap1, Description ==timing[t])
  #   
  #   #if (dim(cor1)[1]==dim(corSamples)[1]) {
  #   #Check if all samples are present in list
  #   indx <- unlist(cor1$Samples )
  #   source=unlist(corSamples$id)
  # 
  #   if  (all(source %in% indx)) {
  #   
  #   t_cor$core=unique(cor1$core)
  #   t_cor$Samples=unique(cor1$Group)
  #   t_cor$Evolution="Common"
  #   tab_cor=rbind(tab_cor,t_cor)
  # 
  #   } else {
  #     
  #     t_cor$core=unique(cor1$core)
  #     t_cor$Samples=unique(cor1$Group)
  #     t_cor$Evolution="unCommon"
  #     tab_cor=rbind(tab_cor,t_cor)
  #   
  #   }
  # }
  #  
  # #Remove row with NA or zero
  # tab_cor=tab_cor[!(is.na(tab_cor$core) | tab_cor$core==""), ]
  # 
  # if ("Common" %in% tab_cor$Evolution) {
  #   tab_cor1=tab_cor[tab_cor$Evolution %in% "Common", ]
  # } else {
  #   tab_cor1=tab_cor
  # }
  # 
  #Select the first 3 columns of cor==my_data
  tab_cor1=cor[,1:3]
  
  cluster=strings_to_BinaryProfiles(df=tab_cor1,CG_pos = CG_pos, class=class)
  cluster1=cluster[,4:(dim(cluster)[2])]
  rownames(cluster1)=1:dim(cluster1)[1]
  cluster1[] <- lapply(cluster1, function(x) as.numeric(as.character(x)))
  cluster1["tot",]=colSums(cluster1)
  cluster1["Freq",]=round(cluster1["tot",]/(dim(cluster1)[1]-1),digits = 2)
   
  # CpGs present in the 80th percentile
  #cores_cutOff1= quantile(x=as.numeric(cluster1["Freq",]))
  clust1=abs(cluster1["Freq",]) >= cores_cutOff[v]
  cores1=paste(colnames(cluster1)[clust1],collapse = "-")
  #Remove empty
  cores1=cores1[cores1 != ""]
  
  # if lenght cores1 =1 reduce percentile di 1
  cores1=unlist(strsplit(cores1, "[-]"))
    
  if (identical(cores1, NULL)) {
    v1=v
    cores1=""
  } else {
    v1=v
    cores1=cores1
  }
  
  
  cores_cutOff1=paste((cores_cutOff[v1]),"%",sep="")
  cores=list(cores=cores1,Evolution=cluster$Evolution,Description=cluster$Description,quantile=cores_cutOff1)
  
  
  return(cores)
  }
  
