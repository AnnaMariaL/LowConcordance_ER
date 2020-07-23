########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: neutral simulations for an empirical p-value distribution of the FET test, empirical data

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


#format.fet: formate the FET result
format.fet<-function(y){
  y$V6<-as.character(y$V6)
  z<-sapply(strsplit(y$V6,"="),function(x) as.numeric(x[2]))
  return(z)
}

#read.fet.data: read FET formats
read.fet.data<-function(x,empirical=F,autosome=T){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes)
  if(!empirical) {
  print("Simulated values!")
  z<-format.fet(y)
  } else {
    print("Empirical values")
    if(autosome){
      y<-subset(y,V1!="X")
      print("Autosomes")
      print(nrow(y))
    } else {
      y<-subset(y,V1=="X")
      print("X")
      print(nrow(y))
    }
    z<-format.fet(y)
  }
  return(z)
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
r<-c(1:10)
my.res<-data.frame(g=integer(),r=integer(),thres.autosome=numeric(),thres.X=numeric())

for(g.iter in g){

  for(r.iter in r){
    #grep files
    emp.data<-list.files(pattern = glob2rx(paste("F0F",g.iter,"SNP_Rep",r.iter,".sync.fet",sep="")))
    sim.autosome<-list.files(pattern = glob2rx(paste("*_autosome_*_rep",r.iter,"F",g.iter,"_simulated.sync.fet",sep = "")))
    sim.X<-list.files(pattern = glob2rx(paste("*_MeanX_*_rep",r.iter,"F",g.iter,"_simulated.sync.fet",sep = "")))
    print(emp.data)
    print(sim.autosome)
    print(sim.X)
  
    #calculate FDR thresholds
    sim.data.autosome<-read.fet.data(sim.autosome,empirical = F,autosome = T)
    sim.data.X<-read.fet.data(sim.X,empirical = F,autosome = F)
    obs.data.autosome<-read.fet.data(emp.data,empirical = T,autosome = T)
    obs.data.X<-read.fet.data(emp.data,empirical = T,autosome = F)
    my.fdr.autosome<-compute_FDR(obs.data.autosome,sim.data.autosome,0.05)
    my.fdr.X<-compute_FDR(obs.data.X,sim.data.X,0.05)
    
    #store results
    temp.res<-data.frame(g=as.integer(g.iter),r=as.integer(r.iter),thres.autosome=my.fdr.autosome,thres.X=my.fdr.X)
    print(paste("Generation:",g.iter,"Replicate",r.iter,":"))
    print(paste("FDR autosome:",my.fdr.autosome,"FDR X:",my.fdr.X))
    my.res<-rbind(my.res,temp.res)
    }
  }

saveRDS(my.res,file = "FDR_FET.rds")
