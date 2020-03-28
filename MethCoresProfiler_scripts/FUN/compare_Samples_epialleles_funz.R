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

compare_Samples_epialleles <- function(LIST,Random_summ1,map) {
  
  em=data.frame(core=NA,Tot=-1)
  em$core=as.character(em$core)

  merge.by.time <- function(a,b) {
    merge(a,b, by=c("core","Tot"),all=TRUE)
  }
  
  em=foreach(h=1:length(LIST), .combine= merge.by.time) %dopar% {
    
    #for (h in 1:length(LIST)) {
    #scomponi il nome del file "i"
    h1=strsplit(LIST[h], "\\.")[[1]] 
    h2=strsplit(h1, "\\_")[[1]]
    
    #if (LIST[h]==Random_summ1) {
     # h3n=5
    #} else {
      h3n=3
    #}
     
      #h3=h2[3]
      h3=paste(h2[6:(length(h2)-h3n)],collapse ="_")
      #carica il file "clones" in base al nome
      summary=read.table(LIST[h], header=TRUE, sep="\t",stringsAsFactors=FALSE)
      #summary=summary[summary$Selection %in% "Positive Selection, pvalue<=0.001",]
      summary1=summary[,c(8,4,5)]
      
      #elimina la r dal nome del campione per fare il merge tra il summary ed il metamap
      sam=h3

      #seleziona un determinato campione dal file metamap
        description=subset(map, select=c(sam))
        description1=t(description)
      
      #Assign gropu name
      colnames(summary1)[1]="core"
      colnames(summary1)[2]="Tot"
      colnames(summary1)[3]=sam
      #clones1$rERR1697396=clones1$rERR1697396*100
      data.frame(summary1)
     #merge(em,summary1,by=c("core","Tot"),all=T)
      #em=merge.by.time(a=em,b=summary1)
     #rm(h)
  }

  em <- em [!is.na(em$core),]
  em[is.na(em)] <- 0
  em[,c(3:(dim(em)[2]))]=em[,c(3:(dim(em)[2]))]
  
  #Conta elementi cores
  #m$dim_core=as.character(str_count(m$core, "-") + 1)
  #Aggiungi collonna Random_control by match with Expected_Freq
  #m <- m %>% merge(Expected_Freq, by="dim_core") %>% select(-dim_core, dim_core = Freq) 
  #colnames(m)[1+(dim(map)[2])]="Random_control"
  #m$Random_control=m$Random_control*
  
  return(em)
}
