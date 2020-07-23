########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: estimate effective population size

#original code: Neda Barghi "Genetic redundancy fuels polygenic adaptation in Drosophila", adapted by Anna Maria Langmüller 

########################
#PACKAGES
########################
rm(list=ls())
#Nest_1.1.8.tar.gz
library(Nest)
########################
#FUNCTIONS
########################

#ne.per.rep: obtain Ne estimates for replicate r at time-point g
ne.per.rep<-function(y,r,g,pop_census,excludechr4=NULL){
  #obtain ancestral allele frequency + coverage
  id.af.0<-which(colnames(y)==paste("F0R",r,sep=""))
  id.cov.0<-which(colnames(y)==paste("F0R",r,"_cov",sep=""))
  #obtain evolved allele frequency + coverage
  id.af.t<-which(colnames(y)==paste("F",g,"R",r,sep=""))
  id.cov.t<-which(colnames(y)==paste("F",g,"R",r,"_cov",sep=""))
  #print(colnames(y)[c(id.af.0,id.af.t,id.cov.0,id.cov.t)])
  
  #estimate Ne
  my.ne<-estimateWndNe(chr=as.character(y$chr), pos=y$pos, wndSize=1000,p0 = y[,id.af.0],pt = y[,id.af.t], cov0=y[,id.cov.0],covt=y[,id.cov.t],t=g,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
  
  #optional filtering
  if(!is.null(excludechr4)) {
    print("Excluding chr. 4")
    my.ne<-subset(my.ne,chr!="4")
    id<-which(my.ne$SNPs<1000 | is.na(my.ne$Np.planI) | my.ne$Np.planI<0)
    if(length(id)>0) {my.ne<-my.ne[-id,]}
  }
  
  my.ne.autosomes<-median(my.ne$Np.planI[my.ne$chr!="X"],na.rm = T)
  id.x<-which(my.ne$chr=="X")
  if(length(id.x)>0){
    my.ne.X<-median(my.ne$Np.planI[my.ne$chr=="X"],na.rm = T)
    my.res<-c(my.ne.autosomes,my.ne.X)
  } else {
    my.res<-my.ne.autosomes
  }
  return(my.res)
}

#my.ne: wrapper function to call ne.per.rep for all replicates
my.ne<-function(x,my.gen,excludechr4=NULL){
  my.data<-read.data(x)
  ne.estimates<-list()
  for(i in c(1:10)) {
    ne.estimates[[i]]<-ne.per.rep(y = my.data,r = i,g = my.gen,pop_census = 300,excludechr4 = excludechr4)
    }
  return(ne.estimates)
}

########################
# ANALYSIS
########################
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/src_sim"
setwd(my_wd)

my.files<-list.files(pattern = glob2rx("*_CovSampled.txt"),path="../data_sim")
my.generations<-seq(10,60,10)

for(i in c(1:length(my.files))){
  Nes<-my.ne(my.files[i],my.generations[i],excludechr4 = TRUE)
  Nes.raw<-do.call("rbind",Nes)
  rownames(Nes.raw)<-paste("F",my.generations[i],"R",c(1:10),sep="")
  colnames(Nes.raw)<-c("autosome","X")
  if(i==1) {
    my.Nes<-Nes.raw
  } else {
    my.Nes<-rbind(my.Nes,Nes.raw)
  }
}

saveRDS(my.Nes,file="NeEstimates.rds")
