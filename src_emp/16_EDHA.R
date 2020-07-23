########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Analyse properties of the different haplotype blocks (EDHA vs. non-EDHA), empirical data


rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)

####################
#PACKAGES
####################

####################
#VARIABLES
####################
input_af<-"Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.AFrising.txt"
cmh_candidates<-"candidates-CMH.rds"
fet_candidates<-"candidates-FET.rds"
the.blocks<-"Dsim_Fl_haplotypeblocks_075_035.sync"
to_exclude<-c("3.17","2.27","3.21","3.48","3.49","3.54")
edha<-c("3.10","3.38","2.18","2.12","3.45","3.40","3.5","3.35","3.25","3.28")
####################
#FUNCTIONS
####################

#read.data.AF: function to read the starting allele frequencies of the SNPs
read.data.AF<-function(x){
  initial<-read.table(x,header=TRUE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=TRUE,sep="\t",colClasses = classes)
  return(y)
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

#merge.candidates: rbind data-frames from CMH and FET tests
merge.candidates<-function(cmh,fet,g){
  colnames(fet)<-c("chr","pos","P","test","gen","rep")
  colnames(cmh)<-c("chr","pos","ref","P","test","gen")
  cmh$rep<-rep(0,nrow(cmh))
  cmh<-subset(cmh,select = c("chr","pos","P","test","gen","rep"))
  print(nrow(cmh))
  print(nrow(fet))
  y<-rbind(cmh,fet)
  y<-subset(y,gen==g)
  return(y)
}

#read.data.blocks: read extracted block information
read.data.blocks<-function(x){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes,na.strings = "-")
  colnames(y)<-c("chr","pos","ID075","ID035")
  return(y)
}

#compute median AF for all blocks (n=99)/replicate and TP
calculateAF<-function(fs,fblock){
  fsM<-merge(fs,fblock,by=c("chr","pos"))
  togrep<-paste("F",rep(seq(0,60,10),each=10),"R",rep(c(1:10),7),sep="")
  fsM<-fsM[,colnames(fsM)%in%togrep]
  return(colMedians(as.matrix(fsM),na.rm = T))
}

s_calculation_approachA<-function(f_s,freqs,thres){
  f_inc_pattern<-increasing(freqs,0,60,thres)
  f_iter<-unique(f_s$blockid)
  f_gens<-c(60,20)
  for (f_gens_iter in f_gens){
    for (f_i in c(1:length(f_iter))){
      f_iter_df<-subset(f_s,blockid==f_iter[f_i])
      f_iter_df<-subset(f_iter_df,gen==f_gens_iter)
      f_iter_inc<-f_inc_pattern[which(rownames(f_inc_pattern)==f_iter[f_i]),]
      f_iter_reps<-which(f_iter_inc==TRUE)
      f_iter_df<-subset(f_iter_df,rep%in%f_iter_reps)
      print(f_iter_df)
      f_iter_res<-data.frame(blockid=unique(f_iter_df$blockid),avgs=mean(f_iter_df$s,na.rm = T),gen=f_gens_iter)
      print(f_iter_res)
      ifelse(f_i==1 & f_gens_iter==60, my_res<-f_iter_res, my_res<-rbind(my_res,f_iter_res))
    }
  }
  return(my_res)
}

#replicate specific increase of blocks
increasing<-function(m,gs,ge,thres){
  find<-paste("F",gs,"R",c(1:10),sep="")
  find<-c(find,paste("F",ge,"R",c(1:10),sep=""))
  mS<-m[,colnames(m)%in%find]
  print(dim(mS))
  mdelta<-mS[,c(11:20)]-mS[,c(1:10)]
  print(head(mdelta))
  mdec<-apply(mdelta, c(1,2), function(x) x>thres)
  return(mdec)
}



####################
#DATA PREP
####################
#ALLELE FREQUENCIES

y<-read.data.AF(input_af)
cmh<-readRDS(cmh_candidates)
fet<-readRDS(fet_candidates)
my_candidates<-merge.candidates(cmh,fet,60)
my_candidates<-subset(my_candidates,select = c("chr","pos"))
my_candidates<-my_candidates[!duplicated(my_candidates),]

#obtain marker SNPs
blocks<-read.data.blocks(the.blocks)
blocks<-blocks[!is.na(blocks$ID035),]
rm(the.blocks)
shared_candidates<-merge(my_candidates,blocks,by=c("chr","pos"))
shared_candidates<-subset(shared_candidates,select = c("chr","pos"))


#create ID and main chromosome coordinates
blocks<-add.main.chr.info(blocks)
blocks$id<-paste(blocks$mainChr,blocks$ID035,sep=".")
length(unique(blocks$id))
blocks$id2<-paste(blocks$chr,blocks$pos,sep=".")

block.lengths<-ddply(blocks,.(id),summarize,min=min(mainPos,na.rm = T),max=max(mainPos,na.rm = T))
block.lengths$length<-block.lengths$max-block.lengths$min
block.lengths$mainChr<-sapply(strsplit(block.lengths$id,"[.]"),"[",1)
block.lengths
#exclude blocks with too low re-discovery rate
block.lengths<-block.lengths[-which(block.lengths$id%in%to_exclude),]

#calculate median allele frequency for each block
blockMedian<-matrix(nrow=nrow(block.lengths),ncol=70)

for(i in c(1:nrow(block.lengths))){
  print(block.lengths$id[i])
  tempBlocks<-subset(blocks,id==block.lengths$id[i])
  #additional line that makes sure the SNPs are also in my candidates
  helper<-merge(tempBlocks,shared_candidates,by=c("chr","pos"))
  tempBlocks<-helper
  rm(helper)
  print(nrow(tempBlocks))
  blockMedian[i,]<-calculateAF(y,tempBlocks)
  print(i)
}
rownames(blockMedian)<-block.lengths$id 
colnames(blockMedian)<-paste("F",rep(seq(0,60,10),each=10),"R",rep(c(1:10),7),sep="")

bdf<-data.frame(blockid=block.lengths$id,min=block.lengths$min,max=block.lengths$max,mainChr=block.lengths$mainChr)
bdf$edha<-F
bdf$edha[bdf$blockid%in%edha]<-T
bdf$p0<-rowMeans(blockMedian[,c(1:10)])

#RECOMBINATION RATE
my.r<-read.table("dsim.rr_LDjump-LOESS-0.1-sex-mimicree2.txt",header = F,skip = 1,sep = "\t")
head(my.r)
my.r$V1<-as.character(my.r$V1)
#change into data frame mainChr, mainPos, r 
#windows are 100 000 bp 
temp<-strsplit(my.r$V1,":")
my.r$chr<-sapply(temp,"[",1)
temp<-sapply(temp,"[",2)
temp2<-strsplit(temp,"[..]")
my.r$pos<-as.numeric(sapply(temp2,"[",3))
my.r<-subset(my.r,select=c("chr","pos","V3"))
my.r<-add.main.chr.info(my.r)
head(my.r)

#for a block: which r-windows does it overlap?
my.recombination.rate<-c()

for(i in c(1:nrow(block.lengths))){
  x<-block.lengths[i,]
  temp.r<-subset(my.r,mainChr==x$mainChr)
  id.s<-which(temp.r$mainPos-100000>=x$min)
  id.e<-which(temp.r$mainPos<=x$max)
  id<-intersect(id.s,id.e)
  if(length(id)==0) id<-c(tail(id.e,1),head(id.s,1))
  my.recombination.rate<-c(my.recombination.rate,mean(temp.r$V3[id]))
}

bdf$r<-my.recombination.rate

#LENGTH
bdf$length<-bdf$max-bdf$min
head(bdf)

#S-ESTIMATE
#estimating selection coefficients for the blocks in single replicates

myne<-readRDS("NeEstimates.rds")
x.blocks<-grep("X",rownames(blockMedian))

s.df<-data.frame(blockid=character(),gen=numeric(),rep=integer(),Ne=numeric(),s=numeric(),p.value=numeric(),p0=numeric())

#for each REP
for(rep.iter in c(1:10)){
  subFreq<-blockMedian[,paste("F",seq(0,60,10),"R",rep.iter,sep="")]
  #for each TP
  for(t.iter in seq(10,60,10)){
    print(paste("R",rep.iter,"F",t.iter,"...."))
    #for each block
    for(block.iter in c(1:nrow(blockMedian))){
      #determine Ne
      useNe<-myne[paste("F",t.iter,"R",rep.iter,sep=""),]
      if(block.iter%in%x.blocks){
        useNe<-useNe[2]
      } else {
        useNe<-useNe[1]
      }
      useNe<-round(useNe)
      names(useNe)<-NULL
      #estimate s + no approximation + correction for census trajectory + p.value 
      my.s<-estimateSH(subFreq[block.iter,paste("F",seq(0,t.iter,10),"R",rep.iter,sep="")],Ne=useNe,haploid = FALSE,h = 0.5,t = seq(0,t.iter,10),N.ctraj = 1000,simulate.p.value = T)
      #create temporary data frame
      temp.df<-data.frame(blockid=rownames(blockMedian)[block.iter],gen=t.iter,rep=rep.iter,Ne=useNe,s=my.s$s,p.value=my.s$p.value,p0=my.s$p0)
      #save in global data frame
      s.df<-rbind(s.df,temp.df)
    }
  }
}

avg_s_approachA_010<-s_calculation_approachA(s.df,blockMedian,0.10)
avg_s_approachA_015<-s_calculation_approachA(s.df,blockMedian,0.15)
avg_s_approachA_020<-s_calculation_approachA(s.df,blockMedian,0.20)
avg_s_approachA_005<-s_calculation_approachA(s.df,blockMedian,0.05)
avg_s_approachA<-rbind(avg_s_approachA_005,avg_s_approachA_010,avg_s_approachA_015,avg_s_approachA_020)
avg_s_approachA$edha<-F
avg_s_approachA$edha[avg_s_approachA$blockid%in%edha]<-T
avg_s_approachA$thres<-rep(c(0.05,0.10,0.15,0.20),each=186)

#calculate selection coefficient ratio 
x<-subset(avg_s_approachA,gen==20)
x<-subset(x,select = c("blockid","avgs","edha","thres"))
y<-subset(avg_s_approachA,gen==60)
y<-subset(y,select = c("blockid","avgs","edha","thres"))
z<-merge(x,y,by=c("blockid","edha","thres"))
head(z)
z$sr<-z$avgs.x/z$avgs.y
colnames(z)<-c("blockid","edha","thres","s_F20","s_F60","sr")

bdf<-merge(bdf,z,by=c("blockid","edha"),all.x = T,all.y = T)

#RISING REPLICATES
for(j in  c(0.05,0.10,0.15,0.20)){
  temp<-increasing(blockMedian,0,20,j)
  temp<-data.frame(blockid=rownames(temp),nrep_20=rowSums(temp))
  temp$thres<-j
  if(j==0.05) temp_res<-temp
  if(j>0.05) temp_res<-rbind(temp_res,temp)
  
}
bdf<-merge(bdf,temp_res,by=c("blockid","thres"))

for(j in  c(0.05,0.10,0.15,0.20)){
  temp<-increasing(blockMedian,0,60,j)
  temp<-data.frame(blockid=rownames(temp),nrep_60=rowSums(temp))
  temp$thres<-j
  if(j==0.05) temp_res<-temp
  if(j>0.05) temp_res<-rbind(temp_res,temp)
  
}

#FINAL DATA SET:
bdf<-merge(bdf,temp_res,by=c("blockid","thres"))

#################
#STATS
#################
#p-value calculations (EDHA features)
sub<-subset(bdf,thres==0.10)
p.af<-wilcox.test(sub$p0~sub$edha)$p.value
p.l<-wilcox.test(sub$length~sub$edha)$p.value
p.r<-wilcox.test(sub$r~sub$edha)$p.value
sub<-subset(bdf,thres==0.05)
sF20.005<-wilcox.test(sub$s_F20~sub$edha)$p.value
sF60.005<-wilcox.test(sub$s_F60~sub$edha)$p.value
sr.005<-wilcox.test(sub$sr~sub$edha)$p.value
reps.005<-wilcox.test(sub$nrep_20~sub$edha)$p.value

sub<-subset(bdf,thres==0.10)
sF20.010<-wilcox.test(sub$s_F20~sub$edha)$p.value
sF60.010<-wilcox.test(sub$s_F60~sub$edha)$p.value
sr.010<-wilcox.test(sub$sr~sub$edha)$p.value
reps.010<-wilcox.test(sub$nrep_20~sub$edha)$p.value

sub<-subset(bdf,thres==0.15)
sF20.015<-wilcox.test(sub$s_F20~sub$edha)$p.value
sF60.015<-wilcox.test(sub$s_F60~sub$edha)$p.value
sr.015<-wilcox.test(sub$sr~sub$edha)$p.value
reps.015<-wilcox.test(sub$nrep_20~sub$edha)$p.value

sub<-subset(bdf,thres==0.20)
sF20.020<-wilcox.test(sub$s_F20~sub$edha)$p.value
sF60.020<-wilcox.test(sub$s_F60~sub$edha)$p.value
sr.020<-wilcox.test(sub$sr~sub$edha)$p.value
reps.020<-wilcox.test(sub$nrep_20~sub$edha)$p.value

my_p<-c(p.af,p.l,p.r,sF20.005,sF20.010,sF20.015,sF20.020,sF60.005,sF60.010,sF60.015,sF60.020,sr.005,sr.010,sr.015,sr.020,reps.005,reps.010,reps.015,reps.020)
names(my_p)<-c("af","length","recombinationrate","sF20-005","sF20-010","sF20-015","sF20-020","sF60-005","sF60-010","sF60-015","sF60-020","sr-005","sr-010","sr-015","sr-020","r-005","r-010","r-015","r-020")
my_p_adj<-p.adjust(my_p,method="BH")
my_p_adj[my_p_adj<0.05]
