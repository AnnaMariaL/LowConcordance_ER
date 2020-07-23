########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: FDR computation, CMH tests, empirical data 

#original code: Neda Barghi "Genetic redundancy fuels polygenic adaptation in Drosophila", adapted by Anna Maria Langmüller 
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)
########################
#PACKAGES
########################
options(stringsAsFactors = FALSE)
library(data.table)

########################
#FUNCTIONS
########################
#read.data: read tab separated data with header
read.data<-function(x,header=TRUE){
  initial<-read.table(x,header=header,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=header,sep="\t",colClasses = classes)
  return(y)
}


#compute_FDR: computes FDR corrected q-values. emp_pval: -log10(pvalue) of empirical data,sim_pval: -log10(pvalue) of simulated data, thresh: FDR threshold, e.g. for 5% provide 0.05
compute_FDR <- function(emp_pval,sim_pval,thresh){
  ecdf_sim <- ecdf(na.omit(sim_pval))
  ecdf_emp <- ecdf(na.omit(emp_pval))
  fdr <- (1-ecdf_sim(emp_pval)) / (1-ecdf_emp(emp_pval))
  fdr[emp_pval > max(sim_pval,na.rm = T)] <- 0
  fdr[sim_pval < min(emp_pval)] <- 1
  fdr[is.na(fdr)]<-1
  candidates <- emp_pval[fdr <= thresh]
  FDR <- min(candidates, na.rm = TRUE)
  return(FDR)
}

########################
#ANALYSIS
########################
g<-seq(60,10,-10)
my.res<-data.frame(g=integer(),thres.autosome=numeric(),thres.X=numeric())
my.ne<-readRDS("NeEstimates.rds")

for(g.iter in g){
  #calculate mean Ne to find the files again
  temp.ne<-my.ne[which(rownames(my.ne)%in%paste("F",g.iter,"R",c(1:10),sep="")),]
  temp.ne.autosome<-round(mean(temp.ne[,1]))
  temp.ne.X<-round(mean(temp.ne[,2]))
  
  sim.file.autosome<-list.files(pattern = glob2rx(paste("*Ne",temp.ne.autosome,"_autosome*simulated.sync.cmh",sep="")))
  sim.file.X<-list.files(pattern = glob2rx(paste("*Ne",temp.ne.X,"_X*simulated.sync.cmh",sep="")))
  emp.file<-list.files(pattern=glob2rx(paste("*F0F",g.iter,"SNP.sync.cmh",sep="")))
  
  print(sim.file.autosome)
  print(sim.file.X)
  print(emp.file)

  sim.data.autosome<-read.data(sim.file.autosome,header = FALSE)
  obs.data<-read.data(emp.file,header = FALSE)
  cmhNull<-sim.data.autosome$V24
  cmhObs<-obs.data$V24[obs.data$V1!="X"]
  
  my.fdr.autosome<-compute_FDR(cmhObs,cmhNull,0.05)
  
  sim.data.X<-read.data(sim.file.X,header=FALSE)
  cmhNull<-sim.data.X$V24
  cmhObs<-obs.data$V24[obs.data$V1=="X"]
  my.fdr.X<-compute_FDR(cmhObs,cmhNull,0.05)
  
  print(paste("Generation",g.iter))
  print(paste("Autosome",my.fdr.autosome))
  print(paste("X",my.fdr.X))
  my.res<-rbind(my.res,data.frame(g=as.integer(g.iter),thres.autosome=my.fdr.autosome,thres.X=my.fdr.X))
  rm(cmhNull,cmhObs,sim.data.X,sim.data.autosome,obs.data,sim.file.autosome,sim.file.X,emp.file)
  print("############")
}

saveRDS(my.res,file="FDR_CMH.rds")


