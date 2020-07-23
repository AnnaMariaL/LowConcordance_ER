########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Haplotype Discovery Rate, empirical data
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)

########################
#PACKAGES
########################
options(stringsAsFactors = FALSE)

library(data.table)
library(poolSeq) #v.0.3.2
library(plyr)
library(factoextra)
library(ggplot2)
library(ggfortify)
############################
#FUNCTIONS
############################

#read.data.blocks: read extracted block information
read.data.blocks<-function(x){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes,na.strings = "-")
  colnames(y)<-c("chr","pos","ID075","ID035")
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

#pol.candidates: get rising/falling allele for a set of snps, given a sync file
pol.candidates<-function(f.x,f.g,f.c){
  f.y<-readRDS(f.x)
  a.y<-alleles(f.y)
  f.c.s<-subset(f.c,gen==f.g)
  a.y<-merge(a.y,f.c.s,by=c("chr","pos"))
  return(a.y)
}

#rediscovered: ratio of SNPs that are rediscovered in each block
rediscovered<-function(source,early,late,blockID){
  source<-subset(source,id==blockID)
  SNPsCount<-nrow(source)
  early<-merge(early,source,by=c("chr","pos"))
  early<-subset(early,select=c("chr","pos","rising"))
  early<-early[!duplicated(early),]
  late<-merge(late,source,by=c("chr","pos"))
  late<-subset(late,select = c("chr","pos","rising"))
  late<-late[!duplicated(late),]
  commons<-merge(early,late,by=c("chr","pos"))
  commons<-subset(commons,select = c("chr","pos","rising.x","rising.y"))
  id<-commons$rising.x==commons$rising.y
  return(c(c(sum(id)/nrow(late)),c(nrow(late)/nrow(source)),length(id)))
}



############################
#VARIABLES
############################
system("cut -f 1,2,85,86 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > Dsim_Fl_haplotypeblocks_075_035.sync")
the.blocks<-"Dsim_Fl_haplotypeblocks_075_035.sync"
read.syncs<-F

###########################
#DATA PREP
###########################
blocks<-read.data.blocks(the.blocks)

blocks<-blocks[!is.na(blocks$ID035),]
rm(the.blocks)
#create ID and main chromosome coordinates
blocks<-add.main.chr.info(blocks)
blocks$id<-paste(blocks$mainChr,blocks$ID035,sep=".")
length(unique(blocks$id))

block.lengths<-ddply(blocks,.(id),summarize,min=min(mainPos,na.rm = T),max=max(mainPos,na.rm = T))
block.lengths$length<-block.lengths$max-block.lengths$min
block.lengths$mainChr<-sapply(strsplit(block.lengths$id,"[.]"),"[",1)
block.lengths

#read sync files for single TP and polarize for the rising allele
if(read.syncs){
  my.syncs<-list.files(pattern = glob2rx("F0F*SNP.sync"))
  gens<-seq(10,60,10)
  for(i in c(1:length(gens))){
    print(my.syncs[i])
    print(gens[i])
    print(gsub(".sync","poolSeq032.rising.rds",my.syncs[i]))
    toSave<-gsub(".sync","poolSeq032.rising.rds",my.syncs[i])
    temp.sync<-read.sync(my.syncs[i],gen = rep(c(0,gens[i]),each=10),repl = rep(c(1:10),2),polarization = "rising")
    saveRDS(temp.sync,toSave)
  }
} else {
  my.syncs<-list.files(pattern = glob2rx("*poolSeq032.rising.rds"))
  print(my.syncs)
}

#obtain list of candidates (t=seq(10,60,10))
#fet<-readRDS("20190416-candidates-FET.rds")
#cmh<-readRDS("20190416-candidates-CMH.rds")
fet<-readRDS("candidates-FET.rds")
cmh<-readRDS("candidates-CMH.rds")
cmh$rep<-rep(0,nrow(cmh))
cmh<-subset(cmh,select = c("chr","pos","logCMH","test","gen","rep"))
colnames(cmh)<-c("chr","pos","P","test","gen","rep")
colnames(fet)<-c("chr","pos","P","test","gen","rep")
candidates<-rbind(cmh,fet)
candidates<-arrange(candidates,chr,pos)

###########################
#HADR
###########################

#read sync files for each TP, polarized for the rising allele
candidates.F10<-pol.candidates(f.x = my.syncs[grep("F0F10",my.syncs)],f.g = 10,f.c = candidates)
candidates.F20<-pol.candidates(f.x = my.syncs[grep("F0F20",my.syncs)],f.g = 20,f.c = candidates)
candidates.F30<-pol.candidates(f.x = my.syncs[grep("F0F30",my.syncs)],f.g = 30,f.c = candidates)
candidates.F40<-pol.candidates(f.x = my.syncs[grep("F0F40",my.syncs)],f.g = 40,f.c = candidates)
candidates.F50<-pol.candidates(f.x = my.syncs[grep("F0F50",my.syncs)],f.g = 50, f.c = candidates)
candidates.F60<-pol.candidates(f.x = my.syncs[grep("F0F60",my.syncs)],f.g = 60, f.c = candidates)

rediscoveredF10<-matrix(ncol = 3,nrow = length(block.lengths$id))
rediscoveredF20<-matrix(ncol = 3,nrow = length(block.lengths$id))
rediscoveredF30<-matrix(ncol = 3,nrow = length(block.lengths$id))
rediscoveredF40<-matrix(ncol = 3,nrow = length(block.lengths$id))
rediscoveredF50<-matrix(ncol = 3,nrow = length(block.lengths$id))
rediscoveredF60<-matrix(ncol = 3,nrow = length(block.lengths$id))

for(id.iter in c(1:length(block.lengths$id))){
  rediscoveredF10[id.iter,]<-rediscovered(blocks,candidates.F10,candidates.F60,block.lengths$id[id.iter])
  rediscoveredF20[id.iter,]<-rediscovered(blocks,candidates.F20,candidates.F60,block.lengths$id[id.iter])
  rediscoveredF30[id.iter,]<-rediscovered(blocks,candidates.F30,candidates.F60,block.lengths$id[id.iter])
  rediscoveredF40[id.iter,]<-rediscovered(blocks,candidates.F40,candidates.F60,block.lengths$id[id.iter])
  rediscoveredF50[id.iter,]<-rediscovered(blocks,candidates.F50,candidates.F60,block.lengths$id[id.iter])
  rediscoveredF60[id.iter,]<-rediscovered(blocks,candidates.F60,candidates.F60,block.lengths$id[id.iter])
}

rownames(rediscoveredF10)<-block.lengths$id
rownames(rediscoveredF20)<-block.lengths$id
rownames(rediscoveredF30)<-block.lengths$id
rownames(rediscoveredF40)<-block.lengths$id
rownames(rediscoveredF50)<-block.lengths$id
rownames(rediscoveredF60)<-block.lengths$id

clusterm<-cbind(rediscoveredF10[,1],rediscoveredF20[,1],rediscoveredF30[,1],rediscoveredF40[,1],rediscoveredF50[,1],rediscoveredF60[,1])
colnames(clusterm)<-paste("F",seq(10,60,10),sep="")
saveRDS(clusterm,file = "ratioRediscoveredOverTime.rds")

#removal of clusters with low ratio of re-discovered marker-SNPs (check rediscoveredFx[,2])
id<-which(rownames(clusterm)%in%c("3.17","2.27","3.21","3.48","3.49","3.54"))
clusterm<-clusterm[-id,]

#clustering, HADR
set.seed(123)
hadr<-clusterm[,-6]
fviz_nbclust(hadr, kmeans, nstart = 2,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method")
hadr_kmeans<-kmeans(hadr,5)

table(hadr_kmeans$cluster)
edha<-names(hadr_kmeans$cluster[hadr_kmeans$cluster=="2"])
hadr[edha,]

hadr_pca<-prcomp(hadr,center = T,scale. = T)
hadr<-as.data.frame(hadr)


hadr$edha<-ifelse(rownames(hadr)%in%edha,T,F)

g<-autoplot(hadr_pca,data=hadr,colour='edha')+theme_minimal()

