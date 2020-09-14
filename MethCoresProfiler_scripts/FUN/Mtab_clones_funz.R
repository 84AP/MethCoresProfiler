
Mtab_clones <- function(sam,CG_pos,b4) {
  
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
  tab$clone=apply(tab,1,paste,collapse="")
  tab=cbind(pop=sam, tab)
  
  tab$n <- rowSums(tab[,nomi_cg]==1,na.rm=T) #data1[1,]==2 Booleano     sum(data1[1,]==2)     sum(data1[1,]==2,na.rm=T)
  
  tab_clones=tab[,c("pop","n", "clone")]
  
  #n_ occurrences of each clone in each sample
  clones_pop=as.data.frame(table(tab$pop, tab_clones$clone))
  colnames(clones_pop)=c("pop", "clone", "n_readsCloneForSample")
  
  #n_reads per sample_without removing unmethylated reads
  tab_compact=aggregate(x=tab_clones, by=list(tab_clones$pop), length) 
  #Alternative: tab_compact=as.data.frame(table(as.character(tab2$pop)))
  rownames(tab_compact)=tab_compact$Group.1 #row indexes_gives the row names of tab_compact the column names of Group.1
  clones_pop$tot_reads_sample=tab_compact[ clones_pop$pop, ]$pop #For each pop of clones_pop $ pop I assign the n_reads taken from tab_compact $ pop
  
  #Freq_clone/sample
  clones_pop$FreqCloneForSample=clones_pop$n_readsCloneForSample/clones_pop$tot_reads_sample
  
  #n_meth calculation
  clones_pop$n=str_count(string = clones_pop$clone, pattern="1") # To each clone of clones_pop $ clone I assign the methylation class taken from clones_types $ n_meth
  
  #n_reads_per methylation class per sample
  tab_compact_group=as.data.frame(table(tab_clones$pop, tab_clones$n))
  #Alternative:tab_compact_groupMeth=aggregate(x=tab_clones, by=list(tab_clones$pop, tab_clones$n_meth), length) #ma rimuove gli zero
  colnames(tab_compact_group)=c("pop", "n", "n_readsGroupForSample")
  tab_compact_group$tmp=paste(tab_compact_group$pop, tab_compact_group$n, sep="")
  clones_pop$tmp=paste(clones_pop$pop, clones_pop$n, sep="")
  rownames(tab_compact_group)=tab_compact_group$tmp 
  clones_pop$n_readsGroupForSample=tab_compact_group[ clones_pop$tmp, ]$n_readsGroupForSample  
  clones_pop$tmp=NULL
  
  #cgpos_cloni
  #nomi_cg=sub(pattern = "X", replacement = "", x = nomi_cg)
  clones_pop$id_pos=-1 #initialize a new column of the clones_pop table
  clones_pop$clone=as.character(clones_pop$clone)
  for (i in 1:(dim(clones_pop)[1] )){    #cycle on the lines
    split_genot=as.numeric(substring(clones_pop[i,]$clone, seq(1,nchar(clones_pop[i,]$clone)),seq(1,nchar(clones_pop[i,]$clone))))    
    split_genot=split_genot==1 #I transform split_genot into a boolean vector, where, only at 2 there will be TRUE
    id= paste(nomi_cg[split_genot],collapse="-")
    if (clones_pop[i,]$n==0) {
      clones_pop[i,]$id_pos=0
    }else{
      clones_pop[i,]$id_pos=id
    }    
  }
  return(clones_pop)
} 