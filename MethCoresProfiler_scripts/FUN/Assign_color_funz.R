Assign_color <- function(m1,Tot,pvalue) {
  #Assign colour based on value
  dim2=m1
  matrix=m1
  #dim2=data.frame(m1[,1])
  #rownames(dim2)=as.character(rownames(m1))
  #em1=em
  
  #merge.by.time <- function(a,b) {
  #  merge(a,b, by="row.names",all=TRUE)
  #}
  
  
  #em1=foreach(h=1:length(LIST)) %dopar% {
  #dim2=foreach(rr=1:ncol(m1), .combine= merge.by.time) %dopar% {   
  
  #Progression_bar
  #pb <- txtProgressBar(title = "Progression Bar", min = 0, max = 100, style = 3)
  #dim2=foreach(rr=1:ncol(m1), .combine=cbind) %dopar% {
  #setTxtProgressBar(pb, rr)
  #dim2=foreach(rr=1:ncol(m1)) %dopar% {
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
    
   
    df2=data.frame(cd25neg_ctr_s1=-1,Random=-1,Tot=-1,Exp_chisqua=-1,Exp_value=1)
    df2=df2[0,]
    #clusterExport(nProcessor, df)
    #i=1
    df2=foreach(i=1:nrow(df), .combine=rbind) %dopar% {
   # for (i in 1:nrow(df)) {
      df1=df[i,]
      contigency_table["Selected","observed_Freq"]=(df1[1,1])*df1$Tot
      contigency_table["expected_Freq","observed_Freq"]=(df1$Random)*df1$Tot
      contigency_table["Selected","Tot_Freq"]=Tot-contigency_table["Selected","observed_Freq"]
      contigency_table["expected_Freq","Tot_Freq"]=Tot-contigency_table["expected_Freq","observed_Freq"]
      x=chisq.test(contigency_table)
      df1$Exp_chisqua=x$statistic
      df1$Exp_value=x$p.value
      data.frame(df1)
      #df2=rbind(df2,df1)
    }
    
    df2=df2[match(rownames(matrix), rownames(df2)),]
    
    #clusterExport(nProcessor, df2)
    #Assign selection type
    outDF=foreach(iii=1:nrow(df2), .combine=rbind) %dopar% {
    #mat=foreach(iii=1:nrow(df2), .combine=rbind) %dopar% {
    #for (iii in 1:nrow(df2)) {
      ncrl=df2[iii,]
      if (ncrl$Exp_value<=pvalue & ncrl[,1]>ncrl[,2]) {
        sel=paste("Positive Selection, pvalue<=",pvalue,sep="")
        color="red"
      } else if (ncrl$Exp_value<=pvalue & ncrl[,1]<ncrl[,2]) {
        sel=paste("Negative Selection, pvalue<=",pvalue,sep="")
        color="blue"
      } else {
        sel="No Selection"
        color="white"
      }
      c(rownames(df2)[iii],sel,color)
    }
    
    outDF=as.data.frame(outDF)
    colnames(outDF)=c("epialleles","Selection","colors")
    
    df2$Selection=outDF$Selection[match(rownames(df2), outDF$epialleles)]

    matrix[,rr]=outDF$colors[match(rownames(matrix), outDF$epialleles)]
  }
  
  dim2=as.matrix(matrix)
  #stopCluster(nProcessor) 
  return(dim2)
}