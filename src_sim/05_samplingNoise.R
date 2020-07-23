########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: code to add binomial sampling to the allele frequencies of the simulated data 

########################
#PACKAGES
########################
rm(list=ls())
options(stringsAsFactors = F)
library(plyr)
######################
#FUNCTIONS
######################
#read.data: read tab separated data with header
read.data<-function(x,header=TRUE){
  initial<-read.table(x,header=header,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=header,sep="\t",colClasses = classes)
  return(y)
}


######################
#VARIABLES
######################
my_pattern<-"*AFrising_cov.txt"
my_gens<-seq(10,60,10)

######################
#DATA & ANALYSIS
######################
#coverage
neda_coverage<-read.csv("../data_emp/data_averageCoverage_PLoS.csv",header = T,sep=";",dec = ",")
colnames(neda_coverage)<-c("Generation","Replicate","Coverage")
neda_coverage$Coverage<-round(neda_coverage$Coverage)

my_files<-list.files(pattern = glob2rx(my_pattern),path="../data_sim/",full.names=T)
print(my_files)

#for all TP
for(i in c(1:6)){
  print(my_files[i])
  #read data
  y<-read.data(my_files[i])
  x<-subset
  #store non-sampled AF
  ns_af<-as.matrix(y[,4:23],nrow=nrow(y),ncol=20)
  #Poisson sampling, empirical coverage
  aim_cov<-neda_coverage$Coverage[neda_coverage$Generation%in%c(0,c(i*10))]
  s_cov<-sapply(aim_cov, function(x) rpois(nrow(ns_af),x))
  #Binomial sampling AF
  s_af<-matrix(0,nrow = nrow(y),ncol = 20)
  s_af<-sapply(c(1:20), function(x) rbinom(nrow(s_af),s_cov[,x],ns_af[,x]))
  #create sync file 
  sync_matrix<-matrix("Z",nrow = nrow(s_af),ncol = ncol(s_af))
  sync_matrix<-sapply(c(1:20), function(x) paste(s_af[,x],":0:",s_cov[,x]-s_af[,x],":0:0:0",sep=""))
  write.table(file = gsub("_AFrising_cov.txt","_CovSampled.sync",my_files[i]),x = cbind(y$chr,y$pos,y$base,as.data.frame(sync_matrix)),col.names = F,row.names = F,sep = "\t",quote = F)
  s_af<-s_af/s_cov
  temp_file<-cbind(y$chr,y$pos,y$base,s_af,s_cov)
  colnames(temp_file)<-colnames(y)
  write.table(file = gsub("_AFrising_cov.txt","_AFrising_CovSampled.txt",my_files[i]),x = temp_file,col.names = T,row.names = F,quote = F,sep = "\t")
}