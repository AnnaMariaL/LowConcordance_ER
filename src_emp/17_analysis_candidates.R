########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: similarity calculations on different genomic analysis heriachies, empirical data 
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)
########################
#PACKAGES
########################
options(stringsAsFactors = F)
library(jaccard)
library(dplyr)
library(reshape2)
library(poolSeq)

########################
#FUNCTIONS
########################
#add.main.chr.info: add main chromosome information to a data frame 
length.2L<-21106324
length.3L<-22263419
add.main.chr.info<-function(x,f.length.2L=length.2L,f.length.3L=length.3L){
  x$chr<-as.character(x$chr)
  x$mainChr<-x$chr
  x$mainChr[which(x$chr%in%c("2L","2R"))]<-"2"
  x$mainChr[which(x$chr%in%c("3L","3R"))]<-"3"
  x$mainPos<-x$pos
  x$mainPos[which(x$chr=="2R")]<-x$mainPos[which(x$chr=="2R")]+f.length.2L
  x$mainPos[which(x$chr=="3R")]<-x$mainPos[which(x$chr=="3R")]+f.length.3L
  x$id<-paste(x$chr,x$pos,sep=".")
  return(x)
}

generate_vector<-function(x){
  x<-subset(x,select=c("chr","pos"))
  x<-paste(x$chr,x$pos,sep=".")
  x<-x[!duplicated(x)]
  return(x)
}

process_data<-function(x){
  y<-readRDS(x)
  z<-generate_vector(y)
  return(z)
}

#Jaccard Index calculation 
determine_J_matrix<-function(x){
  y<-unlist(x)
  y<-unique(y)
  z<-matrix(0,ncol=6,nrow=length(y))
  rownames(z)<-y
  for(i in c(1:ncol(z))){
    z[which(rownames(z)%in%x[[i]]),i]<-1
  }
  print(table(rowSums(z)))
  return(z)
}

determine_jaccard<-function(x,y){
  j<-jaccard(x,y)
  return(j)
}

determine_p_jaccard<-function(x,y,b=10000){
  p0<-jaccard.test(x,y,method = "bootstrap",B=b)
  return(p0$pvalue)
}

#determine candidates for the single time-points
determine_candidates<-function(x,g){
  y<-subset(x,gen==g)
  y$id<-paste(y$chr,y$pos,sep=".")
  y$logPval<-as.numeric(y$logPval)
  y_single<-aggregate(logPval~id,y,function(k) k[which.max(abs(k))])
  y<-merge(y,y_single, by=c("id","logPval"))
  head(y)
  y<-arrange(y,chr,pos)
  return(y)
}

determine_windows<-function(x){
  y<-readRDS(x)
  z<-apply(y[,-ncol(y)],1, function(x) quantile(x,0.99))
  w<-mapply(function(x,y) return(x>=y), x=y[,ncol(y)],y=z)
  w<-w[which(w==TRUE)]
  w<-names(w)
  return(w)
}



############################
#ANALYSIS
############################
generate_candidate_files<-T
generate_data<-T
calculate_jaccard<-T

#extract candidates per generation
if(generate_candidate_files){
  #merge CMH and FET candidates
  cmh<-readRDS("candidates-CMH.rds")
  fet<-readRDS("candidates-FET.rds")
  cmh<-subset(cmh,select = c("chr","pos","test","gen","logCMH"))
  colnames(cmh)<-c("chr","pos","test","gen","logPval")
  cmh$rep<-0
  head(fet)
  colnames(fet)<-c("chr","pos","logPval","test","gen","rep")
  fet<-fet[,c("chr","pos","test","gen","rep","logPval")]
  str(fet)
  str(cmh)
  candidates<-rbind(cmh,fet)

  candidates_F10<-determine_candidates(candidates,10)
  candidates_F20<-determine_candidates(candidates,20)
  candidates_F30<-determine_candidates(candidates,30)
  candidates_F40<-determine_candidates(candidates,40)
  candidates_F50<-determine_candidates(candidates,50)
  candidates_F60<-determine_candidates(candidates,60)

  saveRDS(candidates_F10,file = "candidates.F10.rds")
  saveRDS(candidates_F20,file = "candidates.F20.rds")
  saveRDS(candidates_F30,file = "candidates.F30.rds")
  saveRDS(candidates_F40,file = "candidates.F40.rds")
  saveRDS(candidates_F50,file = "candidates.F50.rds")
  saveRDS(candidates_F60,file = "candidates.F60.rds")
}

#generate the unique lists of candidates 
if(generate_data){
  my.files<-list.files(pattern = glob2rx("*candidates.F*.rds"))
  snp_list<-lapply(my.files, process_data)
  str(snp_list)
  saveRDS(snp_list,"SNPs-candidates-list.rds")

  #Windows
  my.files<-list.files(pattern = glob2rx("Enrichment*.rds"))
  my.files
  window_list<-lapply(my.files, determine_windows)
  str(window_list)
  saveRDS(window_list,"WINDOW-candidates-list.rds")
}

#calculate Jaccard index and the pairwise significance
#p-values corrected with Benjamini Hochberg 

if(calculate_jaccard){
#common alleles
  windows<-readRDS("WINDOW-candidates-list.rds")
  JM<-determine_J_matrix(windows)
  rare<-readRDS("FixedSNPs-F0toF60.rds")
  pairs<-combn(c(1:6),2)
  my_J<-apply(pairs,2,function(x) determine_jaccard(JM[,x[1]],JM[,x[2]]))
  my_J_df<-data.frame(V1=pairs[1,],V2=pairs[2,],V3=my_J)
  my_Jp<-apply(pairs,2,function(x) determine_p_jaccard(JM[,x[1]],JM[,x[2]]))
  my_J_df$pvalue<-my_Jp
  my_J_df$pvalueadj<-p.adjust(my_J_df$pvalue,method="BH")
  my_J_df[my_J_df$pvalueadj<0.05,]


  snps<-readRDS("SNPs-candidates-list.rds")
  JM<-determine_J_matrix(snps)
  #restrict to segregating SNPs
  JM<-JM[-which(rownames(JM)%in%rare),]
  my_J<-apply(pairs,2,function(x) determine_jaccard(JM[,x[2]],JM[,x[1]]))
  my_J_df2<-data.frame(V1=pairs[2,],V2=pairs[1,],V3=my_J)
  my_Jp<-apply(pairs,2,function(x) determine_p_jaccard(JM[,x[2]],JM[,x[1]]))

  my_J_df2$pvalue<-my_Jp
  my_J_df2$pvalueadj<-p.adjust(my_J_df2$pvalue,method="BH")
  my_J_df2[my_J_df2$pvalueadj<0.05,]
  my_J_df$test<-"WINDOW"
  my_J_df2$test<-"SNP"
  my_J<-rbind(my_J_df,my_J_df2)

  my_J$V1<-my_J$V1*10
  my_J$V1<-as.factor(my_J$V1)
  my_J$V2<-my_J$V2*10
  my_J$V2<-as.factor(my_J$V2)
  saveRDS(my_J,file = "Jaccard-SNPs-Windows-COMMON.rds")
}