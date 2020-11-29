## MethCoresProfiler
MethCoresProfiler-master

=============================================
PACKAGE NAME ---> MethCoresProfiler-master
VERSION ----> 2.0
LICENCE---> GNU GPLv3
CREATED BY-> Antonio Pezone <antonio.pezone[at]unina.it>, <antoniopezone[at]gmail.com>
==============================================

INDEX:

-> Requirements
-> About this package
-> Installation
-> Files and their Proper Directories
-> Playing this plugin
-> Known Bugs & Issues
-> Version History
-> Incompatibilities & Save game warnings
-> Credits & Usage


==============================================
REQUIREMENTS:
==============================================
Before running MethCoresProfiler-master scripts make sure following R packages are installed:

require("vegan")
require("Hmisc")
require("data.table")
require("FunChisq")
require("psych")
require("PerformanceAnalytics")
require("gtools")
require("ggpubr")
require("gridExtra")
require("grid")
require("rcompanion")
require("plyr")
require("tidyr")
require("dplyr")
require("entropy")
require("tidyverse")
require("plotrix")
require("igraph")
require("network")
require("foreach")
require("doParallel")
require("stringr")
require("stringi")
require("ggcorrplot")
require("corrplot")
require("scales")
require("ggplot2")
require("reshape2")
require("RelValAnalysis")
require ("RFLPtools")
require ("dendextend")
require("ape")
require("stringdist")
library(reshape)
library("rgl")
require("factoextra")
#require("qdap")
require("stats")
require("rowr")
require("mclust")
require("ComplexHeatmap")
require("gplots")
require("caret")
require("png")
require("circlize")
require("cluster")
library("cowplot")
#require("profr")
require("doSNOW")
require("doMC")
library("progress")
library(pastecs)
library("cowplot")
require("purrr")
require("parallelDist")
require("fastcluster")
require("phangorn")

It is possible to install all necessary packages using pacman, docker or the script "Install_Dependencies.R" (recommended).

==============================================
ABOUT THIS PACKAGE:
===============================================
This package contains R scripts and function to run MethCoresProfiler-master analysis.

Three executable scripts are provided with this version of MethCoresProfiler-master:

  - Run_MethCoresProfiler_setup.R:
    -> after setting all the required parameters, it launches all the scripts in succession. 

  - 1_MethCores_Extractor.R:
    -> requires three types of input files: i. a tab delimited text file of EpiHaplotypes in BINARY format; ii. a tab delimited 
       text file containing information on the CpGs position in the sequence (or string), and; iii. a tab delimited text file 
       containing metadata (information) associated with each sample with the following columns: #SampleID, Tissue, Description, 
       Group, Rep, ID. The first two files can be generated by available tools. Generates several tables reporting: i, All CpG 
       methylation profiles (frequency of single mCpG); ii, The tetrachoric correlations of CpGs; iii, The co-occurrence of two 
       mCpGs; iv, The taxonomic distribution of methylated species; v, The Shannon-Entropy index and the summary/sample, statistics 
       and plots.
  - 2_MethCores_Combinator.R:
    -> takes as input a list of regions, all significant 2-mCpGs combinations annotated by MethCores_Extractor.R and
       a series of filtering parameters. The frequency of all mCpGs combinations (3 mCpGs or more) will be compared 
       to the expected frequency for independent events (x-square for independence statistics, p value ≤ 10-9) and 
       will be annotated; iii. The epialleles containing all significant mCpGs combinations will be extracted; iv. 
       The structure and frequency of significant mCpGs combinations (cores) and individual epialleles will be reported 
       (ComposCore and Tab_Epialleles, respectively). In the absence of significant mCpGs combinations, all epialleles 
       will be brought to the next step
  - 3_MethCores_Analysis.R:
    -> takes as input epialleles annotated by MethCores_Combinator.R, a list of regions and a series of filtering parameters. 
       Generates generating heatmaps of their structure and frequency. In this step, the frequency of each epiallele in the 
       sample is compared (x-square test) to a random control (R), generated automatically by Module 1, with the same number 
       of CpGs and depth (theoretical distribution). Individual significant epialleles are marked in a complex heatmap format 
       (red = high frequency, blue = low frequency and white = neutral). Principal Component Analysis (PCA) of significant 
       epialleles is also performed in this step. Finally, for each sample, the MethCores Analyst extracts the most frequent 
       mCpGs shared by significant epialleles using a decreasing percentage scale (frequency of CpGs in only significant epialleles) 
       starting from 0.9 (minimum 2-mCpGs). If significant epialleles were not found at this stage, the epiallele with the highest 
       frequency will be annotated and its structure will be considered as a stable signature or core.
       The MethCores_Analyst generates three types of indices: i. MethCores Index i.e., the frequency of the methylated cores in 
       the population; ii. Clonality Index, i.e., the frequency of methylated cores normalized to the average methylation in the 
       population; iii. in-phase CpGs Index or Stability Index or Entanglement Index, or E, which measures the rate of coupling or 
       association of at least 2 mCpGs (mCpG1 and mCpG2) in the core 

===============================================
INSTALLATION:
===============================================
MethCoresProfiler-master does not require installations but there are a lot of dependencies.


PREREQUISITES:
===============================================
MethCoresProfiler-master depends on R packages to perform all types of analysis and statistics.

Users can install each dependency:
  1) following R instructions.
  2) Install all dependencies using Install_Dependencies.R
 
MethCoresProfiler-master SCRIPTS:
===============================================
Extract all files from the archive, open the Run_MethCoresProfiler_setup file as text, change the PATH of MethCoresProfiler-master by entering the PATH of the MethCoresProfiler-master folder in and execute the Run_MethCoresProfiler_setup script with Rstudio or R.

===============================================
FILES USED AND THEIR PROPER DIRECTORIES:
===============================================
With this distribution are provided five R scripts and thirty functions. The scripts are in the MethCoresProfiler_scripts folder.

===============================================
USING MethCoresProfiler:
===============================================
Before using MethCoresProfiler-master set the path in the Run_MethCoresProfiler_setup.R:

## Add PATH where is located MethCoresProfiler-master folder
MethCoresProfiler=("/path-to-MethCoresProfiler-master/")

## Add name of folder in which you have put "cgpos", "meta.map" and one folder with all tab delimited "text file of EpiHaplotypes in BINARY format (file.out)" 
## Default is "testData/" in MethCoresProfiler-master folder
Data="testData/"

## Add name of folder containing text file of EpiHaplotypes in BINARY format (file.out) 
## Default is "out/" in MethCoresProfiler-master/testData folder
input="out/"

## Use Rstudio to Run_MethCoresProfiler_setup.R

For detailed instructions on the options and input/output formats of the MethCoresProfiler-master scripts and how to change the parameters we remand the user to the "MethCoresProfiler-master_help.doc" documentation contained in this package.

===============================================
KNOWN ISSUES OR BUGS:
===============================================

===============================================
VERSION HISTORY
===============================================
MethCoresProfiler 2.0