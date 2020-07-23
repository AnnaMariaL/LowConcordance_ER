########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Calculate robustness of ranking in candidate windows, empirical data
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)

########################
#PACKAGES
########################
library(dplyr)
options(stringsAsFactors = F)

####################
#FUNCTIONS
####################
#sig.windows: determine significance state of a window based on simulation percentile
sig.windows<-function(m,thres){
  sims<-apply(m[,-ncol(m)],1,function(x) quantile(x,thres))
  r<-m[,ncol(m)]>=sims
  return(r)
}

keep.extremeP<-function(x){
  x$id<-paste(x$chr,x$pos,sep=".")
  x$pval[which(!is.finite(x$pval))]<-0
  #x$pval[which(is.na(x$pval))]<-0
  x<-arrange(x,id,-pval)
  y<-x[!duplicated(x$id),]
  y<-subset(y,select=c("chr","pos","pval"))
  return(y)
}

#function to format CMH test
read_format_cmh<-function(x){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes)
  y<-subset(y,select=c("V1","V2","V24"))
  colnames(y)<-c("chr","pos","pval")
  y$pval[is.na(y$pval)]<-0
  return(y)
}

#function to format FET test
read_format_fet<-function(x){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes)
  y<-subset(y,select=c("V1","V2","V6"))
  colnames(y)<-c("chr","pos","pval")
  y$pval<-as.character(y$pval)
  y$pval<-gsub(pattern = "1:2=","",y$pval)
  y$pval<-as.numeric(y$pval)
  y$pval[is.na(y$pval)]<-0
  return(y)
}

store_extremes<-function(f_g){
  cmh_pattern<-paste("*_F0F",f_g,"SNP.sync.cmh",sep = "")
  cmh_file<-list.files(pattern = glob2rx(cmh_pattern))
  f_cmh<-read_format_cmh(cmh_file)
  print(cmh_file)
  print(head(f_cmh))
  print("Read cmh")
  for(rep in c(1:10)){
    fet_pattern<-paste("*F0F",f_g,"_Rep",rep,".sync.fet",sep="")
    fet_file<-list.files(pattern = glob2rx(fet_pattern))
    print(fet_file)
    f_fet<-read_format_fet(fet_file)
    print(head(fet_file))
    print(paste("Read",fet_file))
    if(rep==1){
      f_temp<-rbind(f_cmh,f_fet)
      my_res<-keep.extremeP(f_temp)
      rm(f_temp)
    } else {
      my_res<-rbind(my_res,f_fet)
      my_res<-keep.extremeP(my_res)
    }
    print(paste("Filtered",fet_pattern,"for exreme p-values"))
  }
  return(my_res)
}


#determine windows that are enriched in two particular TP 
get_interesting_w<-function(m,s,e){
  ms<-m[,c(s,e)]
  ms<-ms[which(rowSums(ms)>1),]
  return(rownames(ms))
}

#calculate the size of overlap for a certain bin size 
return_overlap<-function(x,y,s){
  x<-x[c(1:s),]
  y<-y[c(1:s),]
  r<-length(intersect(x$id,y$id))
  return(r)
  
}

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

#filter.candidates: extract candidates from a window
#uses Win.Info, add.main.chr.info and Win.filter
filter.candidates<-function(cA,cB,window){
  wW<-Win.Info(window)
  wW<-add.main.chr.info(wW)
  dA<-Win.filter(wca = cA,wdf = wW)
  dB<-Win.filter(wca = cB,wdf = wW)
  dA<-arrange(dA,-P,mainPos)
  dB<-arrange(dB,-P,mainPos)
  return(list(dA,dB))
}


#Win.Info: extract information from window (chr, start, end)
Win.Info<-function(w){
  tw<-unlist(strsplit(w,"-"))
  tw2<-strsplit(tw,"[.]")
  twchr<-sapply(tw2, "[",1)
  twpos<-sapply(tw2,"[",2)
  twdf<-data.frame(chr=twchr,pos=as.numeric(twpos))
  twdf$chr<-as.character(twdf$chr)
  return(twdf)
}


#Win.Filter: filter data frame wca for range in wdf 
Win.filter<-function(wdf,wca){
  wca<-wca[which(wca$mainChr==wdf$mainChr[1] & wca$mainPos>=wdf$mainPos[1] & wca$mainPos<=wdf$mainPos[2]),]
  return(wca)
}

########################
#ANALYSIS
########################

my_gens<-c(10,20,30,40,50,60)
my_thres<-0.99
calculate_overlap<-T


#make .rds files with ALL SNPs and their most extreme -log10(p-values)
for(gen_iter in my_gens) {
  temp<-store_extremes(gen_iter)
  saveRDS(temp,file=paste("F0F",gen_iter,"SNP.pvalues.extremes.rds",sep=""))
}


#determine window matrix 
w<-list.files(pattern = glob2rx("Enrichment*.rds"))

w10<-readRDS(w[1])
w20<-readRDS(w[2])
w30<-readRDS(w[3])
w40<-readRDS(w[4])
w50<-readRDS(w[5])
w60<-readRDS(w[6])

w10_s<-sig.windows(w10,my_thres)
w20_s<-sig.windows(w20,my_thres)
w30_s<-sig.windows(w30,my_thres)
w40_s<-sig.windows(w40,my_thres)
w50_s<-sig.windows(w50,my_thres)
w60_s<-sig.windows(w60,my_thres)

wM<-cbind(w10_s,w20_s,w30_s,w40_s,w50_s,w60_s)
wM<-apply(wM,c(1,2),as.integer)
rownames(wM)<-names(w10_s)
colnames(wM)<-paste("F",seq(10,60,10),sep="")
wM<-wM[-which(rowSums(wM)<1),]


#load SNP sets
g10<-readRDS("F0F10SNP.pvalues.extremes.rds")
g20<-readRDS("F0F20SNP.pvalues.extremes.rds")
g30<-readRDS("F0F30SNP.pvalues.extremes.rds")
g40<-readRDS("F0F40SNP.pvalues.extremes.rds")
g50<-readRDS("F0F50SNP.pvalues.extremes.rds")
g60<-readRDS("F0F60SNP.pvalues.extremes.rds")
g<-list(g10,g20,g30,g40,g50,g60)

#load rare SNPs
rares<-readRDS("FixedSNPs-F0toF60.rds")

#rename + add main chromosomes + filter out rare alleles
for (i in c(1:6)) {
  colnames(g[[i]])<-c("chr","pos","P")
  g[[i]]<-add.main.chr.info(g[[i]])
  g[[i]]<-g[[i]][-which(g[[i]]$id%in%rares),]
  print(i)
}

if(calculate_overlap){
  overlap_df<-data.frame(iter=integer(),overlap=integer(),window=character(),TP=integer())
  #for each intermediate TP 
  for(i in c(1:5)){ 
    #extract interesting windows
    my_windows<-get_interesting_w(wM,i,6)
    print(length(my_windows))
    #for each interesting window
    for(j in my_windows){
      #extract the subsets of the candidates
      w_candidates<-filter.candidates(cA=g[[i]],cB=g[[6]],window = j)
      #store the number of candidates
      w_candidates_nr<-unlist(lapply(w_candidates, nrow))
      if(sum(w_candidates_nr)!=10000) print(paste(j,"in TP",i,"does not have all candidates!"))
      w_candidates_overlap<-sapply(c(1:5000), function(k) return_overlap(x=w_candidates[[1]],y=w_candidates[[2]],s=k))
      w_candidates_overlap_df<-data.frame(iter=c(1:5000),overlap=w_candidates_overlap)
      w_candidates_overlap_df$window<-j
      w_candidates_overlap_df$TP<-i
      overlap_df<-rbind(overlap_df,w_candidates_overlap_df)
      print(paste(j,"from",length(my_windows),"done ..."))
    }
    print(paste("TP",i,"done..."))
  }
  saveRDS(overlap_df,"deltaAUC-allSNPs-segregating.rds")
}


