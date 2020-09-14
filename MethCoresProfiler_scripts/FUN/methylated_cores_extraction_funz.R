methylated_cores_extraction <- function(mer,ttcs,tclust) {
  {
  mer = ttcs[[1]] 
  colnames=list("clone","pop","n_readsCloneForSample","tot_reads_sample","freqCloneForSample","n","n_readsGroupForSample","id_pos") 
  #i=2
  
  for (i in 2:length(ttcs)) {
    mer= merge(mer,ttcs[[i]],by = "clone",all = T)
    
    mer$n_readsCloneForSample.x[is.na(mer$n_readsCloneForSample.x)]=0
    mer$n_readsCloneForSample.y[is.na(mer$n_readsCloneForSample.y)]=0
    mer$n_readsCloneForSample.x=mer$n_readsCloneForSample.x + mer$n_readsCloneForSample.y
    
    mer$tot_reads_sample.x=b4
    
    mer$freqCloneForSample.x[is.na(mer$freqCloneForSample.x)]=0
    mer$freqCloneForSample.y[is.na(mer$freqCloneForSample.y)]=0
    mer$freqCloneForSample.x=mer$freqCloneForSample.x + mer$freqCloneForSample.y
    
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
  mer$freqCloneForSample=mer$freqCloneForSample/length(ttcs)
  mer$n_readsGroupForSample=mer$n_readsGroupForSample/length(ttcs)
  
  # calculate shannon-entropy and normalize to number of CG
  #(1/(dim(CG_pos)[2])/-(sum(mer$freqCloneForSample * log2(mer$freqCloneForSample))))
  mer$shannon_entropy1=(1/(dim(CG_pos)[2])*(entropy(mer$freqCloneForSample, unit="log2"))) #1/nCpG*sum(Entropy)-->calculate population entropoy 
  #entropy.Dirichlet estimates the Shannon entropy H of the random variable Y from the corresponding observed counts y by plug-in of Bayesian estimates of the bin frequencies using the Dirichlet-multinomial pseudocount model.
  mer$Renyi_Diversity=(1/(dim(CG_pos)[2])*(renyi(mer$freqCloneForSample,scales = c(0),hill = T)))
  #entropy.empirical(freq, unit="log2")
  Shannon_entr_tab=mer[1,c("pop","shannon_entropy1")]
  Renyi_entr_tab=mer[1,c("pop","Renyi_Diversity")]
  
  # calculate expected_freq vs observed_freq
  mer$expected_freq=(1/2)^dim(CG_pos)[2]
  mer$observed_freq=mer$freqCloneForSample
  
  #Chi-squared test for given probabilities
  contigency_table=data.frame(matrix(NA, nrow = 2, ncol = 2))
  colnames(contigency_table)=c("observed_freq", "Tot_freq")
  rownames(contigency_table)=c("Selected", "expected_freq")
  
  df=mer[,c(2:5,8,11:12)]
  df2=df
  df2$chisqua=-1
  df2$value=-1
  df2$Fisher_test_pvalue=-1
  df2$Fisher_test_interval=-1
  df2$Fisher_test_conf.level=-1
  df2$Fisher_test_Estimate=-1
  df2=df2[0,]
  
  #i=1
  for (i in 1:nrow(df)) {
    df1=df[i,]
    contigency_table["Selected","observed_freq"]=df1$observed_freq*df1$tot_reads_sample
    contigency_table["expected_freq","observed_freq"]=df1$expected_freq*df1$tot_reads_sample
    contigency_table["Selected","Tot_freq"]=df1$tot_reads_sample-contigency_table["Selected","observed_freq"]
    contigency_table["expected_freq","Tot_freq"]=df1$tot_reads_sample-contigency_table["expected_freq","observed_freq"]
    #contigency_table$RowSum=rowSums(contigency_table[,1:2])
    #contigency_table[3,]=colSums(contigency_table[1:2,])
    #tot_observed_freq=((contigency_table[1,1]/contigency_table[1,3])+(contigency_table[2,1]/contigency_table[2,3]))/2
    #contigency_table=contigency_table[-3,-3]
    x=chisq.test(contigency_table,simulate.p.value = TRUE,B=1000)
    w=wilcox.test(df1$expected_freq, df1$observed_freq, paired=TRUE) 
    #f=fisher.test(contigency_table, conf.level = 0.99,simulate.p.value = TRUE,B = 10000,alternative="one.sided")
    #P=poisson.test(contigency_table, T = 1, r = 1,alternative = c("two.sided"),conf.level = 0.95)
    df1$chisqua=x$statistic
    df1$value=x$p.value
    df1$wilcox_test_pvalue=w$p.value
    #df1$Fisher_test_pvalue=f$p.value
    #df1$Fisher_test_interval=paste(f$conf.int,collapse ="-")
    #df1$Fisher_test_conf.level=0.99
    #df1$Fisher_test_Estimate=f$estimate
    df2=rbind(df2,df1)
  }
  mer=cbind(mer,df2[,8:(dim(df2)[2])])
  
  colnames(mer)[13:14] <- c('x-squared', 'p-value')
  
  b5=as.character(b4)
  Shannon_entropy_tab[rep,b5]=mer$shannon_entropy1[1]
  Renyi_entropy_tab[rep,b5]=mer$Renyi_Diversity[1]
  
  write.table(x = mer, file = paste(sam, "_mclones.txt", sep=""), quote = F, sep = "\t", row.names = F, col.names = T)
  
  #Draw binomial distribution
  flips1=unique(mer[,c(6:7)])
  flips1$Freq=flips1$n_readsGroupForSample/mer$tot_reads_sample[1]
  n=as.data.frame(seq(0,dim(CG_pos)[2],by=1))
  colnames(n)="n"
  flips2=merge(flips1,n,by="n",all=T)
  flips2[is.na(flips2)] <- 0
  flips2$Samples=sam
  
  #taxa_tab=rbind(taxa_tab,flips2)
  
  png(paste("Taxa_Methylated", "_", sam, ".png", sep=""), width=5*150,height=5*150, res = 100,pointsize=4)
  p<-ggplot(data=flips2, aes(x=as.factor(n), y=Freq)) +
    geom_bar(stat="identity", color="black", fill="white")+
    ylim(0,1)+
    ggtitle(paste("Taxa_Methylated", "_", sam, sep=""))+
    xlab("Taxa Species of n_meth")
  #xlim(0,dim(CG_pos)[2])
  
  #p + ggtitle(paste("Binomial_Distibution", "_", sam, sep=""))
  #p + theme(axis.title.x = element_text(size = rel(0.5), angle = 00))
  plot(p)
  dev.off()
  
  mer2= Reduce("+",tclust)/length(tclust)
  
  #Save in pdf the heatmap of the matrix with pvalue 0.01
  png(paste("co_Meth_Occurence", "_", sam, ".png", sep=""), width=5*300,height=5*300, res = 300,pointsize=8)
  corrplot(mer2, type="upper", order="original", sig.level = 0.01, insig = "blank", tl.srt = 100, title=paste("co_Meth_Occurence", "_", sam,sep=""),  mar=c(0,0,1,0))
  dev.off()
  
  meth_list=list(mer = (mer), mer2 =(mer2))
  
  }
  
return(meth_list)

}

