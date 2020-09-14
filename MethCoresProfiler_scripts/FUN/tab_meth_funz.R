
tab_meth <- function(sam,CG_pos,b4) {
  {
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
    
  nomi_cg=as.character(CG_pos[1,])
  colnames(tab)=nomi_cg
  
  tabm=data.frame(t(colSums(tab)/dim(tab)[1]))
  colnames(tabm)=nomi_cg
}
  return(tabm)
  } 
