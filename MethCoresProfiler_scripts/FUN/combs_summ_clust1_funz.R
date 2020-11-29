combs_summ_clust1=function(cl2,sam,clones) {
  #Search for any combination and apply summary function
  summ_clust1= foreach(k=1:nrow(cl2), .combine = rbind) %dopar% {
    #create an empty df
    summary=data.frame(samples=NA, core=NA, Tot=-1, nEpialCore=-1, FreqEpialCore=-1, nEpialUnici=-1,FreqEpialCoreUnici=-1,FreqEpialMoreFreq=-1)
    #assign name
    summary$samples=sam
    #calculate numberbof epialleles
    summary$Tot=sum(clones$n_readsCloneForSample)  
    combs=tidyr::unite_(cl2[k,], paste(colnames(cl2[k,]), collapse="_"), colnames(cl2[k,]),sep="-")
    colnames(combs)="combs"
    #eliminate duplicates
    combs=vapply(lapply(strsplit(combs$combs, "-"), unique), paste, character(1L), collapse = "-")
    #create a list with the cpg set for each epiallele (you can also put this out of the loop)
    cpg_sets <- strsplit(clones$id_pos,"-")
    #look for indices of sets that include all cpgs of core cl2 [k,] regardless of order
    core_ind <- which(sapply(cpg_sets,function(x) { all(cl2[k,] %in% x)}))
    #extract the core
    core <- clones[core_ind,]
    
  }
  #delete all lines containing the NAs
  summ_clust1[complete.cases( summ_clust1),]
  #eliminate duplicates
  summ_clust1=unique( summ_clust1)
  return(summ_clust1)
} 
