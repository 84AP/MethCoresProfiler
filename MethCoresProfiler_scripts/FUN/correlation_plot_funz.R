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

correlation_plot <- function(tclust) {
  
{
  mer2= Reduce("+",tclust)/length(tclust)
  
  #fix ratio
  if (identical(tclust,Meth_unMethclust)) {
    ratio="1CpG_Meth_2CpG_unMeth"
  } else if (identical(tclust,unMeth_Methclust)) {
    ratio="1CpG_unMeth_2CpG_Meth"
  } else {
    ratio="1CpG_Meth_2CpG_Meth"
  }
  
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat) 
    p.mat
  }
  
  #M <- cor(mer2)
  p.mat <- cor.mtest(mer2)
  
  #Salve in pdf la heatmap della matrice con pvalue 0.01
  if (saveAs=="pdf") {
    pdf(paste(sam,"_",ratio,"_co_Occurence_", gene,".pdf", sep=""),paper = "a4")
  } else {
  png(paste(sam,"_",ratio,"_co_Occurence_", gene,".png", sep=""), width=5*300,height=5*300, res = 300,pointsize=8)
  }
    corrplot(mer2, type="upper", order="original", sig.level = 0.01, insig = "blank", tl.srt = 100, p.mat = p.mat,
  title =paste(sam,"_",ratio,"_co_Occurence_", gene,sep=""), mar=c(0,0,1,0), tl.col = "black",cex.main = 0.8)
  dev.off()

  
  #Remove upper part of correlation matrix
  mer2[upper.tri(mer2)] <- NA
  
  #tab_df1=cbind(which(!is.na(mer2),arr.ind = TRUE),na.omit(as.vector(mer2)))
  #Seleziona i valori della matrice = o > di "y1"
  mcorr1 <-  which((mer2>0) & (!is.na(mer2)),arr.ind=TRUE) #estrai row and col numbers compresi quelli della diagonale
  zcorr1 <- mer2[(mer2>0) & (!is.na(mer2))] #estrai valori associati alla righe e colonne di mcorr compresi quelli della diagonale
  cb1=cbind(mcorr1,zcorr1)
  #Crea un dataframe con i valori estratti ed i nomi delle righe e colonne
  
  tab_df1 <- data.frame(cbind(unlist(cb1), #this creates the desired data frame
                              rownames(mer2)[mcorr1[,1]],
                              colnames(mer2)[mcorr1[,2]]
  ))
  
  #estrai le combinazioni
  tab_list=tab_df1[, 4:5]
  #Ordina in maniera crescente le righe
  tab_list=t(apply(tab_list, 1, sort))
  tab_list1=as.data.frame(tab_list)
  colnames(tab_list1)=c("From","to")
  tab_list1$id=paste(tab_list1$From,tab_list1$to,sep="_")
  tab_list1$Freq=round((as.numeric(as.character(tab_df1$zcorr1))),digits = 2)
  
  write.table(x = tab_list1, file = paste(sam,"_",ratio,"_co_Occurence_Histogram_", gene,".txt",sep=""), quote=F,sep="\t",row.names=F)
  
  if (saveAs=="pdf") {
    pdf(paste(sam,"_",ratio,"_co_Occurence_Histogram_", gene,".pdf", sep=""),paper = "a4")
  } else {
  png(paste(sam,"_",ratio,"_co_Occurence_Histogram_", gene,".png",sep=""), width = 5*300,height = 5*300,res = 300,pointsize = 8) 
  }
  j1=ggplot(tab_list1, aes(x=id, y=Freq,fill=id)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) +      # Thinner lines
    #geom_text(aes(label=Freq), position=position_dodge(width=0.9),vjust=-0.25)+
    xlab(ratio) +
    ylab("Frequency") +
    ylim(0,1)+
    ggtitle(paste(sam,"_",ratio,"_co_Occurence_Histogram", gene,sep="")) + 
    #scale_fill_continuous(guide = guide_legend()) +
    theme_classic() +
    theme(legend.position="right", legend.direction="vertical")+
    guides(fill=guide_legend(ncol=1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"),plot.title = element_text(size=5, face="bold"))
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5,color="darkred"), legend.position = "right", legend.box = "vertical",plot.title = element_text(size=5, face="bold"))
  
  addSmallLegend <- function(myPlot, pointSize = 5, textSize = 5, spaceLegend = 0.1) {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_text(size = textSize), 
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }
  
  # Apply on original plot
  j1=addSmallLegend(j1)
  plot(j1)
  dev.off()
  
}
  
}
