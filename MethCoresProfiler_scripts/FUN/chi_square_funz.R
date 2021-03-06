chi_square <- function(mer1,agaist1) {
{
  #Chi-squared test for given probabilities
  contigency_table=data.frame(matrix(NA, nrow = 2, ncol = 2))
  colnames(contigency_table)=c("observed_Freq", "Tot_Freq")
  rownames(contigency_table)=c("Selected", "expected_Freq")
  
  df=mer1
  df$chisqua=-1
  df$value=-1
  #df2$Fisher_test_pvalue=-1
  #df2$Fisher_test_interval=-1
  #df2$Fisher_test_conf.level=-1
  #df2$Fisher_test_Estimate=-1
  df2=df[0,]
  
  #i=1
  for (i in 1:nrow(df)) {
    df1=df[i,]
    
    if (agaist1=="expected") {
    contigency_table["Selected","observed_Freq"]=df1$observed_Freq*df1$Tot
    contigency_table["expected_Freq","observed_Freq"]=df1$expected_Freq*df1$Tot
    contigency_table["Selected","Tot_Freq"]=df1$Tot-contigency_table["Selected","observed_Freq"]
    contigency_table["expected_Freq","Tot_Freq"]=df1$Tot-contigency_table["expected_Freq","observed_Freq"]
    } else {
      contigency_table["Selected","observed_Freq"]=df1$observed_Freq*df1$Tot
      contigency_table["expected_Freq","observed_Freq"]=df1$Random_Freq*df1$Tot
      contigency_table["Selected","Tot_Freq"]=df1$Tot-contigency_table["Selected","observed_Freq"]
      contigency_table["expected_Freq","Tot_Freq"]=df1$Tot-contigency_table["expected_Freq","observed_Freq"]
    }
    
    x=chisq.test(contigency_table)
    #w=wilcox.test(df1$expected_freq, df1$observed_freq, paired=TRUE) 
    #f=fisher.test(contigency_table, conf.level = 0.99,simulate.p.value = TRUE,B = 10000,alternative="one.sided")
    #P=poisson.test(contigency_table, T = 1, r = 1,alternative = c("two.sided"),conf.level = 0.95)
    df1$chisqua=x$statistic
    df1$value=(x$p.value)
    #df1$wilcox_test_pvalue=w$p.value
    #df1$Fisher_test_pvalue=f$p.value
    #df1$Fisher_test_interval=paste(f$conf.int,collapse ="-")
    #df1$Fisher_test_conf.level=0.99
    #df1$Fisher_test_Estimate=f$estimate
    colnames(df1)[3] <- "nEpialTot"
    df2=rbind(df2,df1)
  }
}
  
  return(df2)
}