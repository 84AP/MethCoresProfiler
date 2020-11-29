##### MethCoresProfiler-masterScript.R
### enter the path of testData
#workDir=("./testData/")
## Add PATH of MethCoresProfiler-master
## Add PATH of MethCoresProfiler-master
MethCoresProfiler=("path-to-/MethCoresProfiler-master/")

## Add name of folder in which you have put "cgpos", "meta.map" and one folder with all tab delimited "text file of EpiHaplotypes in BINARY format (file.out)" 
Data="testData/" # folder or ""
## Add name of folder in which you will save cgpos nand metamap files
metaCGpos="metaCGpos/"
## Add name of folder containing text file of EpiHaplotypes in BINARY format (file.out) 
input="out/"
#scriptsDir="MethCoresProfiler_scripts/"

input_outputDir=paste(MethCoresProfiler,Data,sep="")
scriptsDir=paste(MethCoresProfiler,"MethCoresProfiler_scripts/",sep="")
  
#Parameters
nProcessor1=8
#Bootstrapping number to do
B4="default" #"default" or "fix a number, example:10000"
y1=100
#Set pvalue
pvalue1=0.0000000001
#Set average reads or smallest n reads
Reads1="average" ## "smallest"
Remove1="no"  #yes
min1=300
#Set max dim CG_pos to combination
maxCG_pos1=8
#Set max size of RAM (phylo analysis generates very large file
RAM=350
#Plot All samples togheter ("TRUE" or "FALSE")
Plot_All1=TRUE #FALSE"  #or "TRUE"
#make a correlation between 1CpG and 2CpG == Meth and unMeth
Meth_unMeth_corr1="NO"#"YES" #or "NO
saveAs1="png"#pdf or png
Modify_mono1=FALSE #"TRUE" or "FALSE" #if TRUE modify monomethyl
splitForGene1=TRUE #TRUE or FALSE #"splitForGene",

#Put TRUE if you want a temporally analysis
overtime1=TRUE #"FALSE"  #or "TRUE"

## if change only timing
type_Tissue1="same" #"same" "different" 

## reports average weight of significant epialleles containing core. otherwise it reports the weight for each significant epiallels
population_weight1= FALSE #Default is “FALSE”

Groups1=All #"All"  #FALSE  or #"All" #create a unique plot for all samples
vertex_cutOFF1=0.01
HighComplexity1= FALSE # TRUE or FALSE #if TRUE create all combination until n CpGs and keep epialleles containing combinations, 
                                       #if FALSE create all combination with 2CpGS and keep epialleles containing combinations

rm(list = ls()[!ls() %in% c("RAM","metaCGpos","HighComplexity1","vertex_cutOFF1","MethCoresProfiler","input_outputDir","scriptsDir","input","Data","nProcessor1","B4","y1","pvalue1","Reads1","Remove1",
                            "min1","maxCG_pos1","Plot_All1","Meth_unMeth_corr1","saveAs1",
                            "Modify_mono1","splitForGene1","overtime1", "type_Tissue1","population_weight1","Groups1")])

#### Run all scripts in sequence
start.time <- Sys.time()

source(paste(scriptsDir,"1_MethCores_Extractor.R",sep=""))
source(paste(scriptsDir,"2_MethCores_Combinator.R",sep=""))
source(paste(scriptsDir,"3_MethCores_Analyst.R",sep=""))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

