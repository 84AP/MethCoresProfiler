MethCoresIndex <- function(timing,t,tdim5,tmap1,g,timing1,t1) {
  
    my_data = subset(tdim5, Group ==timing[t])
    my_data$Rep=gsub(".*_","",my_data$Samples) 
    samples=as.vector(unique(my_data$Group))
    
    if (dim(my_data)[1]==0) {
      my_data2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
      cores1=""
    } else {
    
    for (i in 1:length(cores_cutOff)) {
      #Use minimal_common_structure2==not take common epiallels amon replicate but only common CpGs among all epiallels
    cores=minimal_common_structure2(cor=my_data,tmap1=tmap1,timing=timing,t=t,v=i)
    cores1=unlist(strsplit(cores$cores, "-"))
    
    if (length(cores1)<=1) next
    
    my_data1=aggregate( Freq ~ core,  my_data, mean )
    my_data2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
    my_data2=my_data2[0,]
    
    for (mm in 1:nrow(my_data1)) {
      mydata2=data.frame(From=NA,to=NA,weight=-1,quantile=-1,sample=NA)
      mydata1=my_data1[mm,]
      
      # make a list with the cpg set for each epiallele (you can put this one out of the loop as well)
      cpg_sets1 =unlist(strsplit(mydata1$core,"-"))
      
      #cores1
      if  (all(cores1 %in% cpg_sets1)==TRUE) {
        mydata2$From[1]=paste(cores1,collapse = "-")
        mydata2$to[1]=paste(setdiff(cpg_sets1,cores1),collapse = "-")
        mydata2$weight[1]=mydata1$Freq
        mydata2$quantile[1]=cores$quantile
        mydata2$sample=samples
        #mydata2=mydata2[!(is.na(mydata2$From) | mydata2$From==""), ]
      }  else {
        mydata2$From=paste(cpg_sets1,collapse = "-")
        mydata2$to=""
        mydata2$weight=mydata1$Freq
        mydata2$quantile[1]=""
        mydata2$sample=samples
      }
      #replace blank with cores
      my_data2=rbind(my_data2,mydata2)
     
    }
    
    my_data2=my_data2[!(is.na(my_data2$From) | my_data2$From==""), ]
    
    if (Groups[g]=="All") {
      my_data2$Tissue=timing1[t1]
    } else {
      my_data2$Tissue=timing[t]
    }
    
    #Change cores dimension
    #if (sum(my_data2$to != "")<1) next
    #if (sum(my_data2$to != "")>=1) break;

    if (sum(my_data2$From != "")>=1) break;

    }

    }
    my_data2=list(my_data2=my_data2,cores=cores1)
    
  return(my_data2)
}