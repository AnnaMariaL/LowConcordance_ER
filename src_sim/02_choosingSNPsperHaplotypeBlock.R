########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: pick the SNPs for the forward simulations

########################
#PACKAGES
########################
rm(list=ls())
library(plyr)
library(reshape2)

########################
#FUNCTIONS
########################

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

#calculate frequencies of A , T, C, and G

freq.baseOnly<-function(h,N=300){
  #split haplotypes
  temp<-unlist(strsplit(h,""))
  temp<-temp[-which(temp==" ")]
  tmp<-table(temp)
  if(sum(tmp)!=c(2*N)) {
    print(paste("weird pos:",pos))
    return(NA)
    }
  base_freqs<-numeric(length = 4)
  base_freqs[1]<-sum(temp=="A")
  base_freqs[2]<-sum(temp=="T")
  base_freqs[3]<-sum(temp=="C")
  base_freqs[4]<-sum(temp=="G")
  base_freqs<-base_freqs/sum(base_freqs)
  names(base_freqs)<-c("A","T","C","G")
  return(base_freqs)
}

#read.data.snps: read extracted marker SNP information
read.data.snps<-function(x){
  initial<-read.table(x,header=FALSE,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=FALSE,sep="\t",colClasses = classes,na.strings = "-")
  colnames(y)<-c("chr","pos","base","ID075","ID035")
  return(y)
}

########################
#DATA & ANALYSIS
########################
# obtain the marker SNP information
system("cut -f 1-3,85,86 ../data_emp/Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > ../data_emp/Dsim_Fl_markerSNPs_075_035.pos" )

the.snps<-"../data_emp/Dsim_Fl_markerSNPs_075_035.pos"

markerSNPs<-read.data.snps(the.snps)
markerSNPs<-markerSNPs[!is.na(markerSNPs$ID035),]
markerSNPs<-add.main.chr.info(markerSNPs)
markerSNPs$id<-paste(markerSNPs$chr,markerSNPs$pos,sep=".")
markerSNPs$blockID<-paste(markerSNPs$mainChr,markerSNPs$ID035,sep=".")
length(unique(markerSNPs$blockID))


#read your ancestral allele frequencies + selection coefficients
block_info<-read.table("../data_emp/S_SAF/S_SAF_0.1AFC.txt")
block_info$id<-c(rep("2",32),rep("3",56),rep("X",11))
block_info$id<-paste(block_info$id,c(1:32,1:56,1:11),sep=".")
block_info<-subset(block_info,select = c("id","V3","V4"))
colnames(block_info)<-c("blockID","s","af")

#add information to the marker SNPs 
markerSNPs<-merge(markerSNPs,block_info,by="blockID",all.x = T)
head(markerSNPs)

# obtain the SNP frequencies of the derived allele in the mimicree file 
#Haplotype information 
my.file<-"../data_sim/Dsim_FL_189Ancestral_wholeGenome-N300.modified.tab"
sexInfo<-FALSE

#skip first line if you have files with sexInformation line 
if(sexInfo){
  mh<-read.csv(my.file,skip = 1,sep="\t",stringsAsFactors = F,header=F)
} else {
  mh<-read.csv(my.file,sep="\t",stringsAsFactors = F,header=F)
}

#name columns
colnames(mh)<-c("chr","pos","ref","ad","h")
mh$id<-paste(mh$chr,mh$pos,sep=".")
#focus on the SNPs that are part of a cluster
mh_cluster<-mh[which(mh$id%in%markerSNPs$id),]

#calculate the base frequencies
mh_freq<-sapply(mh_cluster$h,freq.baseOnly,simplify = F)
mh_freq<-do.call("rbind",mh_freq)
rownames(mh_freq)<-NULL
rownames(mh_freq)<-mh_cluster$id
mh_Freq<-melt(mh_freq)
colnames(mh_Freq)<-c("id","ATCG","freq")
head(mh_Freq)
mh_cluster<-subset(mh_cluster,select = c("chr","pos","ad","id"))
mh_Cluster<-merge(mh_cluster,mh_Freq,by="id")

#keep only options where the base is in ancestral or derived
temp<-strsplit(mh_Cluster$ad,"/")
temp_ancestral<-sapply(temp,"[",1)
temp_derived<-sapply(temp,"[",2)
id_ancestral<-mh_Cluster$ATCG==temp_ancestral
id_derived<-mh_Cluster$ATCG==temp_derived
id<-id_ancestral|id_derived
table(id)
mh_Cluster<-mh_Cluster[id,]
rm(temp,temp_ancestral,temp_derived,id_ancestral,id_derived)

candidates<-merge(mh_Cluster,markerSNPs,by=c("id","chr","pos"),all.x = T)
candidates$delta_af<-abs(candidates$freq-candidates$af)

#keep the one with the minimal deviation from the estimated frequency

candidates_s<-candidates[1,]
for(i in unique(candidates$blockID)){
  temp<-subset(candidates,blockID==i)
  temp<-temp[order(temp$delta_af),]
  candidates_s<-rbind(candidates_s,temp[1,])
}
candidates_s<-candidates_s[-1,]
summary(candidates_s$delta_af)

#create SNP set: 
clusters<-unique(candidates_s$blockID)
snp_set<-candidates_s[1,]
set.seed(20191024)
for(i in clusters){
  temp<-subset(candidates_s,blockID==i)
  print(nrow(temp))
  temp_id<-sample(x=c(1:nrow(temp)),size=1)
  snp_set<-rbind(snp_set,temp[temp_id,])
}

snp_set<-snp_set[-1,]
#correct for a/d ... d SHOULD be the selected one! 
snp_set$ad_corr<-"R/R"
for(i in c(1:nrow(snp_set))){
  temp_ad<-snp_set$ad[i]
  temp_ad<-unlist(strsplit(temp_ad,"/"))
  snp_set$ad_corr[i]<-paste(temp_ad[-which(temp_ad==snp_set$ATCG[i])],temp_ad[which(temp_ad==snp_set$ATCG[i])],sep = "/")
  }

snp_set<-subset(snp_set,select = c("chr","pos","ad_corr","s"))
snp_set$h<-0.5
write.table(snp_set,file = "../data_sim/SNPSet.snp",quote = F,col.names = F,row.names = F,sep = "\t")
