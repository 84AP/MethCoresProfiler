
tetrachoric_corr <- function(sam,CG_pos,b4) {
  
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
  
    
    mydata=tab
    
    #replaces the column name in the.out file with the name of the sequence CG
    colnames(mydata)=CG_pos[1,]
    
    #Run tetrachoric correlation matrix
    tetra= tetrachoric(mydata,y=NULL,correct=FALSE,smooth=FALSE,global=TRUE,weight=NULL,na.rm=FALSE,delete=F)
    
    
    #Extract matrix from "tetra"
    tetra1 <- as.matrix(tetra$rho) 
   # cat(dim(tetra1))
    

  return(tetra1)
} 
