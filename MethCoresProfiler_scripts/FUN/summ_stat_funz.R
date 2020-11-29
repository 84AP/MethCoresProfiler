#create core percentage calculation function
summ_stat <- function(summ, cr,cmb)  {
  
  #if the size of the combination is = 0, everything else is equal to NA
  if (dim(cr)[1] == 0) {
    summ$core=paste(cmb,collapse = "-")
    summ$nEpialCore=NA
    summ$FreqEpialCore=NA
    summ$nEpialUnici=NA
    summ$FreqEpialCoreUnici=NA
    summ$FreqEpialMoreFreq=NA
    
  } else {
    #if combination size is> 0, calculate combination percentage
    summ$core=paste(cmb,collapse = "-")
    summ$nEpialCore=sum(cr$n_readsCloneForSample)
    summ$FreqEpialCore=(summ$nEpialCore/summ$Tot)
    summ$nEpialUnici=dim(cr)[1]
    summ$FreqEpialCoreUnici=(summ$nEpialUnici/summ$Tot)
    
    cr=cr[order(cr$n_readsCloneForSample, decreasing=TRUE),]
    summ$FreqEpialMoreFreq=cr[1,5]
  }
  return(summ)
}
