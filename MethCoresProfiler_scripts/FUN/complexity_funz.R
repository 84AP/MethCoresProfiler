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

complexity <- function(Epi_list) {
  Complexity=data.frame(matrix(NA,ncol=1,nrow=length(Epi_list)))
  colnames(Complexity)="Complexity"
  
  for (ll in 1:length(Epi_list)) {
    ll1=strsplit(Epi_list[ll], "\\_")[[1]] 
    Complexity$Complexity[ll]=ll1[3]
  }
  Complexity=unique(as.vector(unique(t(Complexity))))
  #Complexity=c("dimetili","tetrametili","octometili")
  
  y <- c("dimethyl","tetramethyl","hexamethyl","octamethyl","decamethyl")
  y=y[1:length(Complexity)]
               
Complexity=Complexity[order(match(Complexity, y))]
  
  return(Complexity)
}