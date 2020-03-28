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

# Convert alphanumeric strings to binary profiles using CG_pos files
strings_to_BinaryProfiles <- function(df,CG_pos,class) {
  
  df$core=as.character((df$core))
  n2=dim(df)[2]
  ####### crea dataframe del df_cores
  #crea una lista di tutte le combinazioni
  combs=df$core
  #crea una colonna profiles
  df$profiles= "1"
  
   cluster=foreach(k=1:nrow(df), .combine=rbind) %dopar% {
     #for (k in 1:nrow(df)) {
     #elimina tutte le righe dal dataframe summary
     #summary=summary[0,]
    #crea un dataframe per ogni core
    pcluster=df[k,]
    #crea una colonna profiles con tutti i profili ottenuti
    pcluster$profiles = paste(abs(CG_pos %in% unlist(strsplit(combs[k],"-"))),collapse = "")
    #summary=rbind(summary,pcluster)
    summary=pcluster
     }
  
  if (class=="Methylated") {
    cluster$profiles=cluster$profiles
  } else if (class=="unMethylated") {
    cluster$profiles=gsub("1","2",cluster$profiles)
    cluster$profiles=gsub("0","1",cluster$profiles)
    cluster$profiles=gsub("2","0",cluster$profiles)
  }
  
  #library(tidyr)
  n= max(nchar(as.character(cluster$profiles)))
  cluster=tidyr::separate(cluster, profiles, into = paste0("y", 1:n), sep = 1:(n - 1))
  n1=(n2+dim(CG_pos)[2])
  
  colnames(cluster)[(n2+1):n1]=CG_pos[1,]
   
  return(cluster)
}