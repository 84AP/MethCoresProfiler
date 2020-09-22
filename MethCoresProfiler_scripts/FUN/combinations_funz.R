
combinations<- function(cl,tab_Stat,gene_destination_folder) {
  {
  #add numbering to lines
  cl$id <- seq.int(nrow(cl))
  
  # replace line names with ID names
  row.names(cl)=cl$id
  
  # delete "id" column
  cl$id <- NULL
  
  #separate elements of the string
  combos <- data.frame(do.call('rbind', strsplit(as.character(cl$core),'-',fixed=TRUE)),stringsAsFactors = F)
  combos= as.data.frame(lapply(combos,as.numeric))
  
  #sorts all items in the list
  combos = as.data.frame(t(apply(combos, 1, sort)))
  
  cl2=combos
  

  summary1= combs_summary1(cl2 = cl2,sam = sam, clones = clones)
  stopImplicitCluster()
  
  summ_clust1= combs_summ_clust1(cl2 = cl2,sam = sam, clones = clones)
  stopImplicitCluster()
  
  #write all core files in a single "summary" file which will contain all the combinations obtained until the end of the script
  #z = percentile(summary1 = summary1)
  #summary1=summary1[summary1$PercEpialCore> z,]
  
  # calculate shannon-entropy and normalize to number of CG
  summary1$shannon_entropy=(1/(dim(CG_pos)[2])*(entropy.empirical((summary1$FreqEpialCore), unit="log2")))
  #epiallelic_shannon_entropy=apply(summary1,1, function(x) (1/(dim(CG_pos)[2])/-((summary1$PercEpialCore)*log2((summary1$PercEpialCore)))))
  summary1$RenyiEntropy=(1/(dim(CG_pos)[2])*(RenyiEntropy(summary1$FreqEpialCore,p=0.5)))
  
  #Calcolate Freq combos
  summary1$expected_Freq=-1
  
  for (row in 1:nrow(summary1)) {
    bb=summary1[row,]
    CpGs=unlist(strsplit(summary1$core[row],"-"))
    #select a particular sample from the metamap file
    samples=sub("", "", sam)
    Stat1=tab_Stat[tab_Stat$Samples %in% samples, ]
    CpGs1=subset(Stat1, select=c(CpGs))
    CpGs1=CpGs1
    CpGs1=CpGs1 %>% mutate(Prod = Reduce(`*`, .))
    summary1$expected_Freq[row]=CpGs1$Prod
  }
  summary1$observed_Freq=summary1$FreqEpialCore

  ####### chi_square versus expected
  agaist="expected"
  df2=chi_square(mer1=summary1,agaist1=agaist)
  
  summary1=cbind(summary1,df2[,13:(dim(df2)[2])])
  
  colnames(summary1)[13:(dim(summary1)[2])] <- c('Exp_xSquared', 'Exp_pvalue')
  rm(df2,agaist)
  
  ####### chi_square versus Random
  #Calculate n_Methyl
  summary1$n=as.character(str_count(summary1$core, "-") + 1)
  
  summary1$Random_Freq=-1
  ###Add random_Freq
  for (row1 in 1:nrow(summary1)) {
    row2= summary1[row1,]
    row2= merge(row2,Expected_Freq,by.x="n", by.y="dim_core",all.x=T)
    summary1$Random_Freq[row1]= row2$Freq
  }
  
  agaist="Random"
  df2=chi_square(mer1=summary1,agaist1=agaist)

  summary1=cbind(summary1,df2[,17:(dim(df2)[2])])
  
  colnames(summary1)[17:(dim(df2)[2])] <- c('Random_xSquared', 'Random_pvalue')
  rm(df2,agaist)
  
 #Assign Selection Type Respect to Exp_pvalue
  for (iii in 1:nrow(summary1)) {
    ncrl=summary1[iii,]
    if (ncrl$Exp_pvalue<=pvalue & ncrl$observed_Freq>ncrl$expected_Freq) {
      summary1$Selection[iii]=paste("Positive Selection, pvalue<=",pvalue,sep="")
    } else {
      if (ncrl$Exp_pvalue<=pvalue & ncrl$observed_Freq<ncrl$expected_Freq) {
        summary1$Selection[iii]=paste("Negative Selection, pvalue<=",pvalue,sep="")
      } else {
        summary1$Selection[iii]="No Selection"
      }
    }
  }
  
  #Count n CpGs in cores
  n_cores=max(unique(as.character(str_count(summary1$core, "-") + 1)))
  n_cores1=Complexity[Complexity$Complexity %in% n_cores,]
  dim_n_cores1=paste(n_cores1$Names,"methyl",sep="")
  
  summary1= summary1[!duplicated(summary1$core),,drop=FALSE]
  summ_clust1= summ_clust1[!duplicated(summ_clust1),,drop=FALSE]
  
  #write file
  write.table(x=summary1, file=paste(gene_destination_folder,"/","ComposCore_",dim_n_cores1,"_",class,"_",sam,"_",Clust_Type,"_" ,gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
  write.table(x=summ_clust1, file=paste(gene_destination_folder,"/","Tab_Epialleli_",dim_n_cores1,"_",class,"_", sam, "_",Clust_Type,"_",gene,".txt", sep=""), row.names=F, quote = F, sep="\t")
  
  #eliminate combination already called
  summary1=summary1[!(summary1$core %in% test$core),]
  
  test2=as.data.frame((summary1$core))
  test2[] <- lapply(test2, as.character)
  colnames(test2)="core"
  #test2[complete.cases(test2), ]
  
  summary2=list(summary1=summary1,summ_clust1=summ_clust1,test2=test2)
 }
return(summary2)
}
