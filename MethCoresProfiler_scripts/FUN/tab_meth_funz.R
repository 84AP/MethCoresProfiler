
tab_meth <- function(sam,CG_pos,b4,input_outputDir,input) {
  {
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
    
  nomi_cg=as.character(CG_pos[1,])
  colnames(tab)=nomi_cg
  
  tabm=data.frame(t(colSums(tab)/dim(tab)[1]))
  colnames(tabm)=nomi_cg
}
  return(tabm)
  } 
