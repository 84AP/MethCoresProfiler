Assign_color <- function(m1,Tot,pvalue) {
  #Assign colour based on value
  dim2=m1
  matrix=m1
  
  #Progression_bar
  #pb <- txtProgressBar(title = "Progression Bar", min = 0, max = 100, style = 3)
  #dim2=foreach(rr=1:ncol(m1), .combine=cbind) %:% {
  #setTxtProgressBar(pb, rr)
  for (rr in 1:ncol(m1)) {
    #Chi-squared test for given probabilities
    contigency_table=data.frame(matrix(NA, nrow = 2, ncol = 2))
    colnames(contigency_table)=c("observed_Freq", "Tot_Freq")
    rownames(contigency_table)=c("Selected", "expected_Freq")
    
    df=as.data.frame(m1[,rr])
    rownames(df)=rownames(m1)
    df$Random=m1[,"Random"]
    colnames(df)[1]=colnames(m1)[rr]
    df$Tot=Tot
    
    df2=df
    df2$Exp_chisqua=-1
    df2$Exp_value=-1
    df2=df2[0,]
    
    #i=1
    for (i in 1:nrow(df)) {
      df1=df[i,]
      contigency_table["Selected","observed_Freq"]=(df1[1,1])*df1$Tot
      contigency_table["expected_Freq","observed_Freq"]=(df1$Random)*df1$Tot
      contigency_table["Selected","Tot_Freq"]=Tot-contigency_table["Selected","observed_Freq"]
      contigency_table["expected_Freq","Tot_Freq"]=Tot-contigency_table["expected_Freq","observed_Freq"]
      x=chisq.test(contigency_table)
      df1$Exp_chisqua=x$statistic
      df1$Exp_value=x$p.value
      df2=rbind(df2,df1)
    }
    
    df2=df2[match(rownames(matrix), rownames(df2)),]
    
    #Assign selection type
    #mat=foreach(iii=1:nrow(df2), .combine=rbind) %dopar% {
    for (iii in 1:nrow(df2)) {
      ncrl=df2[iii,]
      if (ncrl$Exp_value<=pvalue & ncrl[,1]>ncrl[,2]) {
        df2$Selection[iii]=paste("Positive Selection, pvalue<=",pvalue,sep="")
        matrix[iii,rr]="red"
      } else if (ncrl$Exp_value<=pvalue & ncrl[,1]<ncrl[,2]) {
        df2$Selection[iii]=paste("Negative Selection, pvalue<=",pvalue,sep="")
        matrix[iii,rr]="blue"
      } else {
        df2$Selection[iii]="No Selection"
        matrix[iii,rr]="white"
      }
    }
    
    dim2[,rr]=matrix[,rr]
  }
  
  dim2=as.matrix(dim2)
  #stopCluster(nProcessor) 
  return(dim2)
}