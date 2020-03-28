 # MethCoresProfiler is a R-script that provides a simple method to trace and track 
 # cores shared by epiallele families in complex populations. 
 # Copyright (C) 2020 author: Antonio Pezone 
 # email: antoniopezone@gmail.com; antonio.pezone@unina.it

 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # any later version.

 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.

 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
    #Carica il file.out dalla lista files
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
  
  #n_occorrenze di ciascun clone in ciascun campione
  clones_pop=as.data.frame(table(tab$pop, tab_clones$clone))
  colnames(clones_pop)=c("pop", "clone", "n_readsCloneForSample")
  
  #n_reads per campione_senza rimuovere le reads non metilate
  tab_compact=aggregate(x=tab_clones, by=list(tab_clones$pop), length) 
  #Alternativa: tab_compact=as.data.frame(table(as.character(tab2$pop)))
  rownames(tab_compact)=tab_compact$Group.1 #indici di riga_assegna ai nomi di riga di tab_compact i nomi della colonna di Group.1
  clones_pop$tot_reads_sample=tab_compact[ clones_pop$pop, ]$pop #Ad ogni pop di clones_pop$pop assegno il n_reads preso da tab_compact$pop
  
  #Freq_clone/sample
  clones_pop$FreqCloneForSample=clones_pop$n_readsCloneForSample/clones_pop$tot_reads_sample
  
  #calcolo n_meth
  clones_pop$n=str_count(string = clones_pop$clone, pattern="1") #Ad ogni clone di clones_pop$clone assegno la classe di metilzione presa da clones_types$n_meth
  
  #n_reads_per classe di metilazione per campione
  tab_compact_group=as.data.frame(table(tab_clones$pop, tab_clones$n))
  #Alternativa:tab_compact_groupMeth=aggregate(x=tab_clones, by=list(tab_clones$pop, tab_clones$n_meth), length) #ma rimuove gli zero
  colnames(tab_compact_group)=c("pop", "n", "n_readsGroupForSample")
  tab_compact_group$tmp=paste(tab_compact_group$pop, tab_compact_group$n, sep="")
  clones_pop$tmp=paste(clones_pop$pop, clones_pop$n, sep="")
  rownames(tab_compact_group)=tab_compact_group$tmp #necessario il rownames
  clones_pop$n_readsGroupForSample=tab_compact_group[ clones_pop$tmp, ]$n_readsGroupForSample  
  clones_pop$tmp=NULL
  
  #cgpos_cloni
  #nomi_cg=sub(pattern = "X", replacement = "", x = nomi_cg)
  clones_pop$id_pos=-1 #inizializzo una nuova colonna della tabella clones_pop
  clones_pop$clone=as.character(clones_pop$clone)
  for (i in 1:(dim(clones_pop)[1] )){    #cicla sulle righe
    split_genot=as.numeric(substring(clones_pop[i,]$clone, seq(1,nchar(clones_pop[i,]$clone)),seq(1,nchar(clones_pop[i,]$clone))))    #ciascun genotipo in ciascuna riga i del tab ciao viene splittata in elementi numerici di un vettore: 
    split_genot=split_genot==1 #Trasformo split_genot in un vettore booleano, dove, solo in corrispondenza del 2 ci sarà TRUE
    id= paste(nomi_cg[split_genot],collapse="-")
    if (clones_pop[i,]$n==0) {
      clones_pop[i,]$id_pos=0
    }else{
      clones_pop[i,]$id_pos=id
    }    
  }
  return(clones_pop)
} 