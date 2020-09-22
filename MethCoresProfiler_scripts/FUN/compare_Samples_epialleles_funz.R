compare_Samples_epialleles <- function(LIST,Random_summ1,map,input_folder2,gene) {
  
  em=data.frame(core=NA,Tot=-1)
  em$core=as.character(em$core)

  merge.by.time <- function(a,b) {
    merge(a,b, by=c("core","Tot"),all=TRUE)
  }
  
  em=foreach(h=1:length(LIST), .combine= merge.by.time) %dopar% {
    
    #for (h in 1:length(LIST)) {
    #separate name
    h1=strsplit(LIST[h], "\\.")[[1]] 
    h2=strsplit(h1, "\\_")[[1]]
    
    #if (LIST[h]==Random_summ1) {
     # h3n=5
    #} else {
      h3n=3
    #}
     
      #h3=h2[3]
      h3=paste(h2[5:(length(h2)-h3n)],collapse ="_")
      #load file by name
      summary=read.table(paste(input_folder2,gene,"/",LIST[h],sep=""), header=TRUE, sep="\t",stringsAsFactors=FALSE)
      #summary=summary[summary$Selection %in% "Positive Selection, pvalue<=0.001",]
      summary1=summary[,c(8,4,5)]
      
      sam=h3

      #select a particular sample from the metamap file
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
  
  
  return(em)
}
