
co_Meth_occurence <- function(sam,CG_pos,b4) {
 class="Methylated"
  
 if (outs1[out]=="Random") {
   n=dim(CG_pos)[2]
   Random=sample(2^n,b4,replace=TRUE) - 1
   Random <- Random + 2^(n)
   names(Random) <- paste0(1:length(Random))
   experiment <- Random %>% 
     purrr::map(function(x){ dectobin(x)}) %>% 
     do.call(rbind,.) %>% as.data.frame() %>% select(-1)
   colnames(experiment)=CG_pos[1,]
 } else {
   
   experiment <- read.table(outs1[out], quote="\"", stringsAsFactors = F)
 }
 
 #Load tab di grandezza b4
 tab=experiment[sample(x = dim(experiment)[1], size = b4, replace = T),]
 tab[] <- lapply(tab[], as.numeric)
 rownames(tab)=1:dim(tab)[1]
 
    #load file out 
    mydata=tab
    
    #replaces the column name in the.out file with the name of the sequence CG
    colnames(mydata)=CG_pos[1,]
    
    #Construct a matrix of co-occurrences of 1
    coM_matrix=matrix(nrow=dim(CG_pos)[2],ncol=dim(CG_pos)[2])
    colnames(coM_matrix)=CG_pos[1,]
    rownames(coM_matrix)=CG_pos[1,]
    #Calculate co-occurrences
    #q=1
    #q1=1
    
    for (q in 1:ncol(mydata)) {
      mydata1=data.frame(mydata[,q])
      colnames(mydata1)=colnames(mydata)[q]
      for (q1 in 1:ncol(mydata)) {
        mydata2=data.frame(mydata[,q1])
        colnames(mydata2)=colnames(mydata)[q1]
        m_mydata=cbind(mydata1,mydata2)
        m_mydata1 = m_mydata[!(apply(m_mydata, 1, function(y) any(y == 0))),]
        cm_mydata1=dim(m_mydata1)[1]/dim(m_mydata)[1]
        coM_matrix[q,q1]=cm_mydata1
      }
    }
  
  
  #write.table(x = coM_matrix, file = paste("Table_coMeth_",class,"_",gene_name1,"_",sam,".txt",sep=""), quote=F,sep="\t",row.names=F)
  
  return(coM_matrix)
} 
