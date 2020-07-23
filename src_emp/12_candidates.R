########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Determine candidate SNPs (CMH + FET ), empirical data
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)

########################
#PACKAGES
########################
options(stringsAsFactors = FALSE)
library(data.table)
library(stringi)

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

#extract.candidates: extract chr, pos, ref, and CMH-score from the output of popoolation
extract.candidates<-function(x){
  x<-subset(x,candidate==TRUE)
  x<-subset(x,select=c(paste("V",c(1:3,24),sep="")))
  colnames(x)<-c("chr","pos","ref","logCMH")
  return(x)
}

#extract.candidates.FET: extract candidates from data.frame
extract.candidates.FET<-function(x){
  x<-subset(x,candidate==TRUE)
  x<-subset(x,select=c("CHROM","POS","P"))
  return(x)
}



#determine.candidates: obtain CMH candidates with exceed a value fdr.thres
determine.candidates<-function(x,fdr.thres){
  x$candidate<-rep(FALSE,nrow(x))
  x$candidate[which(x$V24>=fdr.thres)]<-TRUE
  return(x)
}

#get.replicate: obtain replicate number from the FET names
get.replicate<-function(x){
  my.r<-stri_extract(x,regex = "Rep[0-9]+")
  my.r<-gsub("Rep","",my.r)
  my.r<-as.integer(my.r)
  return(my.r)
}

########################
#candidate SNPs, CMH-test
########################
fdr.cmh<-readRDS("FDR_CMH.rds")
colnames(fdr.cmh)<-c("gen","autosome","X")
fdr.cmh$thres<-round(fdr.cmh$autosome)
fdr.cmh$thresX<-round(fdr.cmh$X)
fdr.cmh

#empirical CMH tests:
my.cmh<-list.files(pattern = glob2rx("*.cmh"))
my.cmh<-my.cmh[-grep("simulated",my.cmh)]
my.cmh<-rev(my.cmh)
print(my.cmh)

my.candidates.cmh<-data.frame(chr=character(),pos=integer(),ref=character(),logCMH=numeric())
my.generations<-seq(60,10,-10)

for(cmh.iter in c(1:length(my.cmh))){
  print(cmh.iter)
  my.gen<-my.generations[cmh.iter]
  temp.data<-read.data(my.cmh[cmh.iter],header = FALSE)
  temp.data<-determine.candidates(temp.data,fdr.thres = fdr.cmh$thres[which(fdr.cmh$gen==my.gen)])
  temp.candidates<-extract.candidates(temp.data)
  temp.candidates$test<-rep("CMH",nrow(temp.candidates))
  temp.candidates$gen<-rep(my.gen,nrow(temp.candidates))
  print(table(temp.candidates$chr))
  my.candidates.cmh<-rbind(my.candidates.cmh,temp.candidates)
}


saveRDS(my.candidates.cmh,"candidates-CMH.rds")

########################
#candidate SNPs, FET
########################

my.res<-data.frame(g=integer(),r=integer(),thres.autosome=numeric(),thres.X=numeric())


my.fet<-list.files(pattern = glob2rx("*.fet"))
my.fet<-my.fet[-grep("simulated",my.fet)]
my.fet<-rev(my.fet)
print(my.fet)

fet.cutoffs<-readRDS("FDR_FET.rds")

my.candidates.fet<-data.frame(CHROM=character(),POS=integer(),P=numeric(),test=character(),gen=integer(),rep=integer())
my.generations<-rep(seq(60,10,-10),each=10)

for(fet.iter in c(1:length(my.fet))){
  temp.fet<-my.fet[fet.iter]
  my.gen<-my.generations[fet.iter]
  my.rep<-get.replicate(temp.fet)
  temp.cutoff<-fet.cutoffs$thres.autosome[which(fet.cutoffs$g==my.gen & fet.cutoffs$r==my.rep)]
  print(paste("File:",temp.fet))
  print(paste("Replicate:",my.rep))
  print(paste("Generation:",my.gen))
  print(paste("FDR cutoff:",temp.cutoff))
  temp.data<-read.data.FET(temp.fet)
  print(head(temp.data))
  temp.data$candidate<-F
  temp.data$candidate[temp.data$P>=temp.cutoff]<-T
  temp.candidates<-extract.candidates.FET(temp.data)
  temp.candidates$test<-rep("FET",nrow(temp.candidates))
  temp.candidates$gen<-rep(my.gen,nrow(temp.candidates))
  temp.candidates$rep<-rep(my.rep,nrow(temp.candidates))
  print(table(temp.candidates$CHROM))
  my.candidates.fet<-rbind(my.candidates.fet,temp.candidates)
}


saveRDS(my.candidates.fet,"candidates-FET.rds")


