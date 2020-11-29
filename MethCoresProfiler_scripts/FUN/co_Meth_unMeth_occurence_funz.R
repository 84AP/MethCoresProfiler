
co_Meth_unMeth_occurence <- function(sam,CG_pos,b4,input_outputDir,input) {
 #class="Methylated"
  
  if (outs1[out]=="Random") {
    n=dim(CG_pos)[2]
    if (n>=51) {
      experiment=as.data.frame(matrix(ncol=dim(CG_pos)[2],nrow=b4))
      experiment[is.na(experiment)] <- 0
      experiment=map_df(experiment, function(x) {x[sample(c(TRUE, NA), prob = c(0.8, 0.2), size = length(x), replace = TRUE)]})
      # A tibble: 10 x 3
      colnames(experiment)=CG_pos[1,]
      experiment[is.na(experiment)] <- 1
      experiment=as.data.frame(experiment)
      #for (nr in 1:nrow(experiment)) {
      # experiment[nr,nr]=1
      #}
      #experiment[is.na(experiment)] <- 0
    } else {
    Random=sample(2^n,b4,replace=TRUE) - 1
    Random <- Random + 2^(n)
    names(Random) <- paste0(1:length(Random))
    experiment <- Random %>% 
      purrr::map(function(x){ dectobin(x)}) %>% 
      do.call(rbind,.) %>% as.data.frame() %>% select(-1)
    colnames(experiment)=CG_pos[1,]
  }
  } else {
    
    experiment <- read.table(paste(input_outputDir,"/",input,"/",outs1[out],sep=""), quote="\"", stringsAsFactors = F)
  }
  
  #Load tab di grandezza b4
  tab=experiment[sample(x = dim(experiment)[1], size = b4, replace = T),]
  tab[] <- lapply(tab[], as.numeric)
  rownames(tab)=1:dim(tab)[1]
  
   #set value1 and value2
   value1=1
   value2=0
    
   #load
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
        m_mydata1=m_mydata[m_mydata[,1] ==value1, ] 
        m_mydata1=m_mydata1[m_mydata1[,2] ==value2, ]
        #m_mydata1 = m_mydata[!(apply(m_mydata, 1, function(y) any(y == 0))),]
        cm_mydata1=dim(m_mydata1)[1]/dim(m_mydata)[1]
        coM_matrix[q,q1]=cm_mydata1
      }
    }
  
  
  #write.table(x = coM_matrix, file = paste("Table_coMeth_",class,"_",gene_name1,"_",sam,".txt",sep=""), quote=F,sep="\t",row.names=F)
  
  return(coM_matrix)
  
 
  
}
