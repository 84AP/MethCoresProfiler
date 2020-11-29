# Convert alphanumeric strings to binary profiles using CG_pos files
strings_to_BinaryProfiles <- function(df,CG_pos,class) {
  
  df$core=as.character((df$core))
  n2=dim(df)[2]
  ####### creates dataframe of df_cores
   # make a list of all combinations
  combs=df$core
  #create a profiles column
  df$profiles= "1"
  
   cluster=foreach(k=1:nrow(df), .combine=rbind) %dopar% {
     #for (k in 1:nrow(df)) {
     
     #summary=summary[0,]
    #create a df for each core
    pcluster=df[k,]
    #creates a profiles column with all the obtained profiles
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