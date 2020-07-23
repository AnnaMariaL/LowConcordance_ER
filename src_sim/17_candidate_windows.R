########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Determine candidate windows

rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_sim/"
setwd(my_wd)
########################
#PACKAGES 
########################
options(stringsAsFactors = FALSE)
library(data.table)
library(stringi)

#####################
#FUNCTIONS
#####################
#read.data.AF: function to read the starting allele frequencies of the SNPs
read.data.AF<-function(x){
  initial<-read.table(x,header=TRUE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=TRUE,sep="\t",colClasses = classes)
  return(y)
}

#obtain.fixed.SNPs: obtain ID (chr[.]pos) of SNPs that fix/get lost in at least one replicate in a certain generation
obtain.fixed.SNPs<-function(x,g){
  y<-read.data.AF(x)
  y.id<-paste(y$chr,y$pos,sep=".")
  y<-y[,colnames(y)%in%paste("F",g,"R",c(1:10),sep="")]
  print(str(y))
  id.low<-apply(y,1,function(x) any(x==0))
  id.up<-apply(y,1,function(x) any(x==1))
  id<-id.low|id.up
  return(y.id[which(id==TRUE)])
}

#candidate.vector: obtain character vector of candidates from a certain generation
candidate.vector<-function(df,g){
  df<-subset(df,gen==g)
  df<-subset(df,select=c("chr","pos"))
  df<-paste(df$chr,df$pos,sep=".")
  df<-df[!duplicated(df)]
  return(df)
}

#call.sims: call simulation to sample same amount of single bins from whole genome
call.sims<-function(fs=s,fafcanr=af.ca.nr,start,end){
  id<-sapply(c(1:length(fafcanr)), function(i) sample.SNPs(x=fs,j = fafcanr[i]),simplify = T)
  id<-unlist(id)
  ca.sim<-rep(0,length(fs))
  names(ca.sim)<-names(fs)
  ca.sim[names(ca.sim)%in%names(id)]<-1
  nr.ca.sim<-mapply(function(x,y,z=ca.sim) sum(ca.sim[x:y]), x=start,y=end)
  return(nr.ca.sim)
}


#####################
#ANALYSIS
#####################
window<-5000
my.thres<-0.99
epochs<-1000
my.gen<-c(10,20,30,40,50,60)
conservative<-T

#obtain list of candidates (t=seq(10,60,10))
fet<-readRDS("candidates-FET.rds")
cmh<-readRDS("candidates-CMH.rds")
cmh$rep<-rep(0,nrow(cmh))
cmh<-subset(cmh,select = c("chr","pos","logCMH","test","gen","rep"))
colnames(cmh)<-c("chr","pos","P","test","gen","rep")
colnames(fet)<-c("chr","pos","P","test","gen","rep")
candidates<-rbind(cmh,fet)
rm(cmh,fet)

#filter SNP-set for segregating alleles in all TP and REPS --> put boolean conservative to TRUE
if(conservative){
  my.files<-list.files(pattern = "Sampled.txt")
  fixedF0<-obtain.fixed.SNPs(my.files[1],0)
  fixedF10<-obtain.fixed.SNPs(my.files[1],10)
  fixedF20<-obtain.fixed.SNPs(my.files[2],20)
  fixedF30<-obtain.fixed.SNPs(my.files[3],30)
  fixedF40<-obtain.fixed.SNPs(my.files[4],40)
  fixedF50<-obtain.fixed.SNPs(my.files[5],50)
  fixedF60<-obtain.fixed.SNPs(my.files[6],60)
  fixed<-c(fixedF0,fixedF10,fixedF20,fixedF30,fixedF40,fixedF50,fixedF60)
  fixed<-fixed[!duplicated(fixed)]
  saveRDS(fixed,file = "FixedSNPs-F0toF60.rds")
  #remove SNPs that are fixed in any replicate at any generation
  myAF<-read.data.AF(my.files[1])
  id.myAF<-paste(myAF$chr,myAF$pos,sep=".")
  id.myAF.filtered<-id.myAF[-which(id.myAF%in%fixed)]
  
  
  #filter candidates for these SNPs
  candidates$id<-paste(candidates$chr,candidates$pos,sep=".")
  candidates<-subset(candidates,id%in%id.myAF.filtered)
  
  #filter AF from candidates
  myAF<-myAF[-which(id.myAF%in%fixed),]
  myAF<-myAF[,colnames(myAF)%in%paste("F0R",c(1:10),sep="")]
  meanAF<-rowMeans(myAF)
  names(meanAF)<-id.myAF.filtered
  head(meanAF)
  
} else {
  #obtain starting allele frequencies 
  myAF<-read.data.AF("F0F60SNP_AFrising_cov.txt")
  myAF<-subset(myAF,select = c("chr","pos","base",paste("F0R",c(1:10),sep="")))
  meanAF<-rowMeans(myAF[,-c(1:3)])
  names(meanAF)<-paste(myAF$chr,myAF$pos,sep=".")
  }

rm(myAF)


for(g in my.gen){
  
  #obtain candidates 
  tocheck<-candidate.vector(candidates,g)

  

  ###########################
  #ENRICHMENT
  ###########################

  #vector s: numeric vector storing the indices
  #vector ca: numeric vector: 1... candidate at pos in s, 0 ... no candidate at pos in s
  #af .... allele frequencies of the snps
  #afbin .... binned allele frequencies
  #window ... window size for the analysis 
  
  s<-c(1:length(meanAF))
  names(s)<-names(meanAF)

  ca<-rep(0,length(s))
  names(ca)<-names(s)
  ca[names(ca)%in%tocheck]<-1

  af<-meanAF

  afbin<-cut(af,breaks = quantile(af,probs = seq(0,1,0.1)),include.lowest = T)
  names(afbin)<-names(s)

  #create start and end position of each window based on the window size
  #do not overlap between main chromosomes
  id2<-grep("^2[LR][.]",names(s))
  id3<-grep("^3[LR][.]",names(s))
  id4<-grep("^4[.]",names(s))
  idX<-grep("X[.]",names(s))

  a<-numeric()
  b<-numeric()
  la<-seq(id2[1],id2[length(id2)],by = window)
  lb<-la+c(window-1)
  lb[length(lb)]<-max(id2)

  a<-c(a,la)
  b<-c(b,lb)

  la<-seq(id3[1],id3[length(id3)],by = window)
  lb<-la+c(window-1)
  lb[length(lb)]<-max(id3)
  a<-c(a,la)
  b<-c(b,lb)

  la<-seq(id4[1],id4[length(id4)],by = window)
  lb<-la+c(window-1)
  lb[length(lb)]<-max(id4)
  a<-c(a,la)
  b<-c(b,lb)
  
  la<-seq(idX[1],idX[length(idX)],by = window)
  lb<-la+c(window-1)
  lb[length(lb)]<-max(idX)

  a<-c(a,la)
  b<-c(b,lb)

  print(str(a))
  print(str(b))

  #calculate #candidates for each window
  #nr.ca ... number of candidates for each window
  nr.ca<-mapply(function(x,y,z=ca) sum(ca[x:y]),x=a,y=b)

  #obtain afbins of candidates
  af.ca<-afbin[ca!=0]

  #count SNPs in each class
  af.ca.nr<-table(af.ca)
  
  #actual simulations
  system.time(simData<-replicate(epochs,call.sims(start = a,end = b)))

  #name matrix and store the raw matrix, last column:observed counts
  my.names<-paste(names(s)[a],names(s)[b],sep="-")
  rownames(simData)<-my.names
  simData<-cbind(simData,nr.ca)
  tosave<-paste("Enrichment-",epochs,"Epochs-",window,"SNPs-generation",g,".rds",sep="")
  saveRDS(simData,file = tosave)

  #determine if window is enriched or not
  cutoffs<-apply(simData[,-ncol(simData)],1,function(x) quantile(x,my.thres))
  my.res<-simData[,ncol(simData)]>cutoffs
  names(my.res)<-my.names
  print(paste("Generation",g,":"))
  print(table(my.res))
}


