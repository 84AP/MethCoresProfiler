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