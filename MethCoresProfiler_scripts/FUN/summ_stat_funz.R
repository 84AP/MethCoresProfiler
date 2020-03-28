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

#crea funzione calcolo percentuali core
summ_stat <- function(summ, cr,cmb)  {
  
  #se la dimensione della combinazione è = 0, tutto il resto è uguale a NA
  if (dim(cr)[1] == 0) {
    summ$core=paste(cmb,collapse = "-")
    summ$nEpialCore=NA
    summ$FreqEpialCore=NA
    summ$nEpialUnici=NA
    summ$FreqEpialCoreUnici=NA
    summ$FreqEpialMoreFreq=NA
    
  } else {
    #se la dimensione della combinazione è >0 , calcola percentuale combinazione
    summ$core=paste(cmb,collapse = "-")
    summ$nEpialCore=sum(cr$n_readsCloneForSample) #unici e non
    summ$FreqEpialCore=(summ$nEpialCore/summ$Tot)
    summ$nEpialUnici=dim(cr)[1]
    summ$FreqEpialCoreUnici=(summ$nEpialUnici/summ$Tot)
    
    cr=cr[order(cr$n_readsCloneForSample, decreasing=TRUE),]
    summ$FreqEpialMoreFreq=cr[1,5]
  }
  return(summ)
}
