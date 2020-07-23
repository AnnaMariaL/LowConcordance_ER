########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: check the amount of SNPs that have fixed throughout the experiment

rm(list=ls()) 
options(stringsAsFactors = F)
library(poolSeq)
library(plyr)
library(ggplot2)
wd<-"/Volumes/Temp1/early_response_Dsim_FL/data/"
setwd(wd)

#--------------------------
read_sync<-T
if(read_sync){
mysync<-read.sync("Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.forAF.sync",gen = rep(seq(0,60,10),each=10),repl = rep(c(1:10),7),polarization = "rising")
saveRDS(mysync,file="Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.forAF_poolSeq032_rising.rds")
} else {
  mysync<-readRDS("Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.forAF_poolSeq032_rising.rds")
}

#x ... path to candidates
#g ... generation of candidates
#thres .... min/max AF , e.g. thres = 0.01, minAF = 0.01, maxAF=1-0.01
#s ... sync file, all gens and reps
NREPS_FIXED<-function(x,g,thres,s){
  res_df<-data.frame(nreps=integer(),thres=numeric(),cand_gen=integer(),gen=integer(),nsnp=integer()) #empty data.frame 
  y<-readRDS(x) #read data 
  y<-subset(y,select = c("chr","pos","rising"))
  y<-y[!duplicated(y),]
  sa<-alleles(sync = s,chr = y$chr,pos = y$pos) #check allelic state + correct (although not neccessary in this implementation)
  id_polar<-which(sa$rising!=y$rising)
  j<-seq(g+10,60,10) #generations to check 
  
  for(j_iter in j){
    saf<-af(sync = s,chr = y$chr,pos = y$pos,repl = c(1:10),gen =j_iter)
    saf[id_polar,]<-1-saf[id_polar,]
    nreps<-apply(saf,1,function(x) sum(x<=thres | x>=(1-thres)))
    temp_df<-data.frame(nreps=nreps)
    temp_df$thres<-thres
    temp_df$cand_gen<-g
    temp_df$gen<-j_iter
    temp_df$nsnp<-nrow(y)
    res_df<-rbind(res_df,temp_df)
  }
  return(res_df)
}

gen10<-NREPS_FIXED("candidates.F10.rds",10,thres = 0.01,s = mysync)
gen20<-NREPS_FIXED("candidates.F20.rds",20,thres = 0.01,s = mysync)
gen30<-NREPS_FIXED("candidates.F30.rds",30,thres = 0.01,s = mysync)
gen40<-NREPS_FIXED("candidates.F40.rds",40,thres=0.01,s=mysync)
gen50<-NREPS_FIXED("candidates.F50.rds",50,thres=0.01,s=mysync)

gen<-rbind(gen10,gen20,gen30,gen40,gen50)
gen_counts<-aggregate(list(gen$nreps),list(gen$cand_gen,gen$gen,gen$nsnp),function(x) sum(x>7))
head(gen_counts)
colnames(gen_counts)<-c("candgen","gen","ncand","SNP")
gen_counts$SNPr<-gen_counts$SNP/gen_counts$ncand
gen_counts<-arrange(gen_counts,candgen,gen)
gen_counts$label<-round(gen_counts$SNPr,3)
gen_counts

g<-ggplot(data=gen_counts,aes(x=gen,y=SNPr,shape=as.factor(candgen)))+geom_point(size=3)+theme_classic()
g<-g+xlab("generation")+ylab("ratio of fixed SNPs")+theme(text = element_text(size = 15),legend.position = "top")+labs(shape="candidates")
g

