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

 combs_summ_clust1=function(cl2,sam,clones) {
  #Cerca ogni combinazione e applica funzione summary
  summ_clust1= foreach(k=1:nrow(cl2), .combine = rbind) %dopar% {
    #crea un dataframe vuoto
    summary=data.frame(samples=NA, core=NA, Tot=-1, nEpialCore=-1, FreqEpialCore=-1, nEpialUnici=-1,FreqEpialCoreUnici=-1,FreqEpialMoreFreq=-1)
    #assegna nome al campione nel summary
    summary$samples=sam
    #calcola il numero degli epialleli
    summary$Tot=sum(clones$n_readsCloneForSample)  #unici e non
    combs=tidyr::unite_(cl2[k,], paste(colnames(cl2[k,]), collapse="_"), colnames(cl2[k,]),sep="-")
    colnames(combs)="combs"
    #elimina duplicati dalla lista
    combs=vapply(lapply(strsplit(combs$combs, "-"), unique), paste, character(1L), collapse = "-")
    #crea una lista con il set di cpg per ogni epiallele (la puoi mettere pure fuori al loop questa)
    cpg_sets <- strsplit(clones$id_pos,"-")
    #cerca gli indici dei set che includono tutte le cpg del core cl2[k,] a prescindere dall ordine
    core_ind <- which(sapply(cpg_sets,function(x) { all(cl2[k,] %in% x)}))
    #estrai il core
    core <- clones[core_ind,]
    
  }
  #elimina tutte le righe che contengono gli NA
  summ_clust1[complete.cases( summ_clust1),]
  #elimina duplicati
  summ_clust1=unique( summ_clust1)
  return(summ_clust1)
} 
