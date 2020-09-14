cores_extraction <- function(mer,ttcs,tclust) {
  {
    
  mer = ttcs[[1]] 
  
  #set population name
    class="Methylated"
    n_type="n_Meth"
  
  colnames=list("clone","pop","n_readsCloneForSample","tot_reads_sample","FreqCloneForSample","n","n_readsGroupForSample","id_pos") 
  #i=2
  
  for (i in 2:length(ttcs)) {
    mer= merge(mer,ttcs[[i]],by = "clone",all = T)
    
    mer$n_readsCloneForSample.x[is.na(mer$n_readsCloneForSample.x)]=0
    mer$n_readsCloneForSample.y[is.na(mer$n_readsCloneForSample.y)]=0
    mer$n_readsCloneForSample.x=mer$n_readsCloneForSample.x + mer$n_readsCloneForSample.y
    
    mer$tot_reads_sample.x=b4
    
    mer$FreqCloneForSample.x[is.na(mer$FreqCloneForSample.x)]=0
    mer$FreqCloneForSample.y[is.na(mer$FreqCloneForSample.y)]=0
    mer$FreqCloneForSample.x=mer$FreqCloneForSample.x + mer$FreqCloneForSample.y
    
    mer$n_readsGroupForSample.x[is.na(mer$n_readsGroupForSample.x)]=0
    mer$n_readsGroupForSample.y[is.na(mer$n_readsGroupForSample.y)]=0
    mer$n_readsGroupForSample.x=mer$n_readsGroupForSample.x + mer$n_readsGroupForSample.y
    
    mer$id_pos.x[is.na(mer$pop.x)] = mer$id_pos.y[is.na(mer$id_pos.x)]
    #mer$id_pos.x[is.na(mer$id_pos.x)] = mer$id_pos.y[is.na(mer$id_pos.x)]
    mer$n.x[is.na(mer$pop.x)] = mer$n.y[is.na(mer$n.x)]
    
    mer$pop.x=sam
    
    mer=mer[,1:8]
    colnames(mer)=colnames
  }
  #Divide by the number of arrays of tcs
  mer$n_readsCloneForSample=mer$n_readsCloneForSample/length(ttcs)
  mer$FreqCloneForSample=mer$FreqCloneForSample/length(ttcs)
  mer$n_readsGroupForSample=mer$n_readsGroupForSample/length(ttcs)
  # calculate shannon-entropy and normalize to number of CG
  #(1/(dim(CG_pos)[2])/-(sum(mer$FreqCloneForSample * log2(mer$FreqCloneForSample))))
  mer$shannon_entropy1=(1/(dim(CG_pos)[2])*(entropy(mer$FreqCloneForSample, unit = "log2")))
  #entropy.Dirichlet estimates the Shannon entropy H of the random variable Y from the corresponding observed counts y by plug-in of Bayesian estimates of the bin Frequencies using the Dirichlet-multinomial pseudocount model.
  mer$RenyiEntropy=(1/(dim(CG_pos)[2])*(RenyiEntropy(mer$FreqCloneForSample,0.5)))
  #entropy.empirical(Freq, unit="log2")
  Shannon_entropy_tab=mer[1,c("pop","shannon_entropy1")]
  Renyi_entropy_tab=mer[1,c("pop","RenyiEntropy")]
  
  epiallelic_tab=mer[,c(1,5,8)]
  colnames(epiallelic_tab)[2]=unique(mer$pop)
    
  # calculate expected_Freq vs observed_Freq
  #mer$expected_Freq=(1/2)^dim(CG_pos)[2]
  #mer$observed_Freq=mer$FreqCloneForSample
  
  #chi_suqare funz
  #df2=chi_square(mer1=mer)
  
  #df2=df$df2
  #mer=cbind(mer,df2[,8:(dim(df2)[2])])
  
  #colnames(mer)[13:14] <- c('x-squared', 'p-value')
  
  write.table(x = mer, file = paste("./",destination_folder,"/",sam,"_",class,"_",gene,"_clones.txt", sep=""), quote = F, sep = "\t", row.names = F, col.names = T)
  
  #Correlation matrix building
  cor=mer[,c("id_pos","FreqCloneForSample")]
  colnames(cor)[2]=sam
  
  #Draw binomial distribution
  flips1=unique(mer[,c(6:7)])
  flips1$Freq=flips1$n_readsGroupForSample/mer$tot_reads_sample[1]
  n=as.data.frame(seq(0,dim(CG_pos)[2],by=1))
  colnames(n)="n"
  flips2=merge(flips1,n,by="n",all=T)
  flips2[is.na(flips2)] <- 0
  flips2$Samples=sam
  
  taxa_tab=flips2
  taxa_tab_mean =aggregate(taxa_tab$Freq, list(taxa_tab$Samples,taxa_tab$n), mean)
  colnames(taxa_tab_mean)=c("Samples","n","Freq")
  taxa_tab_sd =aggregate(taxa_tab$Freq, list(taxa_tab$Samples,taxa_tab$n), sd)
  colnames(taxa_tab_sd)=c("Samples","n","sd")
  taxa_tab1=merge(taxa_tab_mean,taxa_tab_sd,by=c("Samples","n"),all = TRUE)
  taxa_tab1[is.na(taxa_tab1)] <- 0
  
  if (saveAs=="pdf") {
    pdf(paste("./",destination_folder,"/",sam,"_",class,"_Taxa_",gene,".pdf", sep=""),paper = "a4")
  } else {
  png(paste("./",destination_folder,"/",sam,"_",class,"_Taxa_",gene,".png", sep=""), width=5*150,height=5*150, res = 100,pointsize=4)
  }
    p<-ggplot(data=taxa_tab1, aes(x=as.factor(n), y=Freq)) +
    geom_bar(stat="identity", color="black", fill="white")+
    geom_errorbar(aes(ymin=Freq-sd, ymax=Freq+sd), width=.2,
                  position=position_dodge(.9))+
    ylim(0,1)+
    ggtitle(paste(sam,"_",class,"_Taxa_",gene, sep=""))+
    xlab(paste("Taxa Species of ", n_type,sep=""))
  #xlim(0,dim(CG_pos)[2])
  
  #p + ggtitle(paste("Binomial_Distibution", "_", sam, sep=""))
  #p + theme(axis.title.x = element_text(size = rel(0.5), angle = 00))
  plot(p)
  dev.off()
  
  mer2= Reduce("+",tclust)/length(tclust)
  
  #Save in pdf the heatmap of the matrix with pvalue 0.01
  if (saveAs=="pdf") {
    pdf(paste("./",destination_folder,"/",sam,"_co_",class,"_Occurence_", gene,".pdf", sep=""),paper = "a4")
  } else {
  png(paste("./",destination_folder,"/",sam,"_co_",class,"_Occurence_", gene,".png", sep=""), width=5*300,height=5*300, res = 300,pointsize=8)
  }
    corrplot(mer2, type="upper", order="original", sig.level = 0.01, insig = "blank", tl.srt = 100, title=paste(sam,"_co_",class,"_Occurence_", gene,sep=""),  mar=c(0,0,1,0))
  dev.off()
  
  meth_list=list(mer = (mer), mer2 =(mer2),
  Shannon_entropy_tab=(Shannon_entropy_tab),
  Renyi_entropy_tab=(Renyi_entropy_tab),
  cor=(cor),
  taxa_tab=(taxa_tab1),
  epiallelic_tab=epiallelic_tab)
  
  }
  
return(meth_list)
}

