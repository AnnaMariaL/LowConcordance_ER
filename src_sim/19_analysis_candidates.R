########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: similarity calculations on different genomic analysis heriachies 
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_sim/"
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

###################################
#NR OF CANDIDATES OVER TP
###################################
targets<-read.table("SNPset.snp",skip = 1,sep = "\t",header = F)

colnames(targets)<-c("chr","pos","ad","s","h")
targets<-add.main.chr.info(targets)

candidates.F10<-readRDS("candidates.F10.rds")
candidates.F10<-add.main.chr.info(candidates.F10)
found.F10<-merge(candidates.F10,targets,by=c("mainChr","mainPos","chr","pos","id"))
nrow(found.F10)

candidates.F20<-readRDS("candidates.F20.rds")
candidates.F20<-add.main.chr.info(candidates.F20)
found.F20<-merge(candidates.F20,targets,by=c("mainChr","mainPos","chr","pos","id"))
nrow(found.F20)

candidates.F30<-readRDS("candidates.F30.rds")
candidates.F30<-add.main.chr.info(candidates.F30)
found.F30<-merge(candidates.F30,targets,by=c("mainChr","mainPos","chr","pos","id"))

candidates.F40<-readRDS("candidates.F40.rds")
candidates.F40<-add.main.chr.info(candidates.F40)
found.F40<-merge(candidates.F40,targets,by=c("mainChr","mainPos","chr","pos","id"))

candidates.F50<-readRDS("candidates.F50.rds")
candidates.F50<-add.main.chr.info(candidates.F50)
found.F50<-merge(candidates.F50,targets,by=c("mainChr","mainPos","chr","pos","id"))

candidates.F60<-readRDS("candidates.F60.rds")
candidates.F60<-add.main.chr.info(candidates.F60)
found.F60<-merge(candidates.F60,targets,by=c("mainChr","mainPos","chr","pos","id"))



###############################
#RR : generate Data Frame 

###############################
rr<-read.csv("../data_emp/dsim.rr_LDjump-LOESS-0.1-sex-mimicree2.txt",skip = 1,header = F,sep="\t")

temp<-strsplit(rr$V1,"[..]")
temp2<-sapply(temp,"[",1)
temp3<-strsplit(temp2,":")
rr$chr<-sapply(temp3,"[",1)
rr$start<-sapply(temp3,"[",2)
rr$start<-as.numeric(rr$start)
rr$end<-sapply(temp,"[",3)
rr$end<-as.numeric(rr$end)

head(rr)
rm(temp3,temp2,temp)
rr<-subset(rr,select = c("V3","chr","start","end"))
colnames(rr)<-c("cM","chr","start","end")
head(rr)

#add rr to targets 
targets$rr<-numeric(length = nrow(targets))
for(i in c(1:nrow(targets))){
  temp<-targets[i,]
  temp_rr<-subset(rr,chr==temp$chr)
  id_up<-which(temp_rr$start<=temp$pos)
  id_down<-which(temp_rr$end>=temp$pos)
  id<-intersect(id_up,id_down)
  targets$rr[i]<-temp_rr$cM[id]
}

View(targets)

###########################################
#PROPERTIES OF EARLY VS. NON-EARLY TARGETS
###########################################
#read with poolSeq + polarized for the rising allele 
my.sync<-readRDS("DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x_rising.rds")

early<-subset(found.F20,select = c("chr","pos"))
early<-inner_join(early,subset(found.F30,select = c("chr","pos")),by=c("chr","pos"))
early<-inner_join(early,subset(found.F40,select = c("chr","pos")),by = c("chr","pos"))
early<-inner_join(early,subset(found.F50,select = c("chr","pos")),by=c("chr","pos"))
early<-inner_join(early,subset(found.F60,select = c("chr","pos")),by=c("chr","pos"))
late<-subset(found.F60,select = c("chr","pos"))
late<-anti_join(late,early)

#A:nr.rising replicates
thres<-c(0.05,0.1,0.15,0.2)
gens<-c(20)
for(j in gens){
for(i in thres){
  #number of replicates rising in early  targets
  earlyF0<-af(sync = my.sync,chr = early$chr,pos = early$pos,repl = c(1:10),gen = 0)
  earlyFt<-af(sync = my.sync,chr = early$chr,pos = early$pos,repl = c(1:10),gen = j)
  earlyD<-earlyFt-earlyF0
  early_reps<-apply(earlyD,1,function(x) sum(x>i))
  #number of replicates rising in late targets
  lateF0<-af(sync = my.sync,chr = late$chr,pos = late$pos,repl = c(1:10),gen = 0)
  lateFt<-af(sync = my.sync,chr = late$chr,pos = late$pos,repl = c(1:10),gen = j)
  lateD<-lateFt-lateF0
  late_reps<-apply(lateD,1,function(x) sum(x>i))
  temp_rising<-data.frame(chr=c(early$chr,late$chr),pos=c(early$pos,late$pos),nrep=c(early_reps,late_reps))
  temp_rising$thres<-i
  temp_rising$gen<-j
  temp_rising$type<-c(rep("early",nrow(early)),rep("late",nrow(late)))
  if(j==20 & i==0.05) {
    nrepDF<-temp_rising
  } else {
    nrepDF<-rbind(nrepDF,temp_rising)
  }
}
  }
nrepDF$thres<-as.factor(nrepDF$thres)


#B: s-estimates 
earlyS<-inner_join(targets,early,by=c("chr","pos"))
earlyS<-subset(earlyS,select = c("chr","pos","s"))
lateS<-inner_join(targets,late,by=c("chr","pos"))
lateS<-subset(lateS,select = c("chr","pos","s"))
SDF<-data.frame(chr=c(earlyS$chr,lateS$chr),pos=c(earlyS$pos,lateS$pos),s=c(earlyS$s,lateS$s))
SDF$type<-c(rep("early",nrow(earlyS)),rep("late",nrow(lateS)))

#C: starting frequency
saf<-af(sync = my.sync,chr = targets$chr,pos = targets$pos,gen = 0,repl = c(1:10))
saf<-apply(saf,1,mean)
ADF<-splitLocusID(names(saf))
ADF$af<-saf

earlyAF<-inner_join(ADF,early,by=c("chr","pos"))
lateAF<-inner_join(ADF,late,by=c("chr","pos"))
ADF<-data.frame(chr=c(earlyAF$chr,lateAF$chr),pos=c(earlyAF$pos,lateAF$pos),af=c(earlyAF$af,lateAF$af))
ADF$type<-c(rep("early",nrow(earlyAF)),rep("late",nrow(lateAF)))
head(ADF)

#D: recombination rate 
earlyRR<-inner_join(targets,early,by=c("chr","pos"))
earlyRR<-subset(earlyRR,select = c("chr","pos","rr"))
earlyRR$type<-"early"
lateRR<-inner_join(targets,late,by=c("chr","pos"))
lateRR<-subset(lateRR,select = c("chr","pos","rr"))
lateRR$type<-"late"
RRDF<-rbind(earlyRR,lateRR)
ggboxplot(data=RRDF,x="type",y="rr")
wilcox.test(rr~type,data=RRDF)

#calculate p.values
#p-value calculations (EDHA features)
s<-subset(nrepDF,thres==0.05 & gen==20)
r005<-wilcox.test(nrep~type,data=s)$p.value
s<-subset(nrepDF,thres==0.10 & gen==20)
r010<-wilcox.test(nrep~type,data = s)$p.value
s<-subset(nrepDF,thres=0.15 & gen==20)
r015<-wilcox.test(nrep~type,data = s)$p.value
s<-subset(nrepDF,thres==0.20 & gen==20)
r020<-wilcox.test(nrep~type,data = s)$p.value
ps<-wilcox.test(s~type,data = SDF)$p.value
paf<-wilcox.test(af~type,data = ADF)$p.value
prr<-wilcox.test(rr~type,data = RRDF)$p.value

my_p<-c(r005,r010,r015,r020,ps,paf,prr)
names(my_p)<-c("r-005","r-010","r-015","r-020","s","af","recombinationrate")
my_p_adj<-p.adjust(my_p,method = "BH")
my_p_adj[my_p_adj<0.05]