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
 
co_Meth_occurence <- function(sam,CG_pos,b4) {
 class="Methylated"
  
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
    
    #Costruisci una matrice di co-occorenze di 1
    coM_matrix=matrix(nrow=dim(CG_pos)[2],ncol=dim(CG_pos)[2])
    colnames(coM_matrix)=CG_pos[1,]
    rownames(coM_matrix)=CG_pos[1,]
    #Calcola le co-occorenze
    #q=1
    #q1=1
    
    for (q in 1:ncol(mydata)) {
      mydata1=data.frame(mydata[,q])
      colnames(mydata1)=colnames(mydata)[q]
      for (q1 in 1:ncol(mydata)) {
        mydata2=data.frame(mydata[,q1])
        colnames(mydata2)=colnames(mydata)[q1]
        m_mydata=cbind(mydata1,mydata2)
        m_mydata1 = m_mydata[!(apply(m_mydata, 1, function(y) any(y == 0))),]
        cm_mydata1=dim(m_mydata1)[1]/dim(m_mydata)[1]
        coM_matrix[q,q1]=cm_mydata1
      }
    }
  
  
  #write.table(x = coM_matrix, file = paste("Table_coMeth_",class,"_",gene_name1,"_",sam,".txt",sep=""), quote=F,sep="\t",row.names=F)
  
  return(coM_matrix)
} 
