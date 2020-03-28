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
    #Carica il file.out dalla lista files
    experiment <- read.table(outs1[out], quote="\"", stringsAsFactors = F)
  }
  
  
  #Load tab di grandezza b4
  tab=experiment[sample(x = dim(experiment)[1], size = b4, replace = T),]
  tab[] <- lapply(tab[], as.numeric)
  rownames(tab)=1:dim(tab)[1]
  
    #Carica file out 
    mydata=tab
    
    #sostituisce il nome della colonna nel file.out con il nome delle CG della sequenza
    colnames(mydata)=CG_pos[1,]
    
    #Esegui tetrachoric correlation matrix
    tetra= tetrachoric(mydata,y=NULL,correct=FALSE,smooth=FALSE,global=TRUE,weight=NULL,na.rm=FALSE,delete=F)
    
    
    #Estrai matrice da "tetra"
    tetra1 <- as.matrix(tetra$rho) 
   # cat(dim(tetra1))
    

  return(tetra1)
} 
