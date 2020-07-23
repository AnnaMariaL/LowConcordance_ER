########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: neutral simulations for an empirical p-value distribution of the CMH test, empirical data

#original code: Neda Barghi "Genetic redundancy fuels polygenic adaptation in Drosophila", adapted by Anna Maria Langmüller 
rm(list=ls())
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/src_emp/"
setwd(my_wd)
########################
#PACKAGES
########################
options(stringsAsFactors = FALSE)
library(data.table)
library(plyr)
########################
#FUNCTIONS
########################
#convert.to.sync: converts the counts of simulated drift to sync format. For simplicity all polymorphic sites are A/C combinations.
convert.to.sync <- function(x,y){paste(x, 0,y-x, 0, 0, 0, sep=":")}

# fast.trajectory: simulate drift trajectories
fast.trajectory <- function(p0, Ne, t, s=0, h=0.5) {
  traj <- matrix(NA, ncol=length(t), nrow=max(length(p0), length(s), length(h)), dimnames=list(c(), paste0("F", t)))
  if(0 %in% t)
    traj[,"F0"] <- p0
  
  g <- 1
  p <- p0
  q <- 1-p0
  wAA <- 1+s
  wAa <- 1+h*s
  waa <- 1
  # simulate allele frequencies across time
  while(g <= max(t)) {
    # apply selection and random drift
    p <- (wAA*p^2 + wAa*p*q) / (wAA*p^2 + wAa*2*p*q + waa*q^2)
    if(!is.na(Ne))
      p <- rbinom(length(p), Ne, p)/Ne
    q <- 1-p
    
    # if necessary then save current allele frequency to results
    if(g %in% t)
      traj[,paste0("F", g)] <- p
    
    g <- g+1
  }
  
  if(nrow(traj) == 1)
    return(as.vector(traj))
  
  return(traj)
}


#calculate.average.ancestral: compute the average of starting allele frequencies in 10 replicates of founder popualtion
calculate.average.ancestral<-function(y){
  togrep<-paste("F0R",c(1:10),sep="")
  id.ancestral<-which(colnames(y)%in%togrep)
  print(colnames(y)[id.ancestral])
  y$safDt<-apply(y[,id.ancestral],1,sum)
  y$safDt<-y$safDt/10
  return(y)
}

#calculate.counts: add noise onto allele counts
calculate.counts<-function(y,autosome,my.Ne,my.gen,Bcov,Ecov){
  if(autosome){
    simTraj<-fast.trajectory(p0=y$safDt[y$chr!="X"],Ne=my.Ne*2,t=c(0,my.gen))
  } else {
    simTraj<-fast.trajectory(p0=y$safDt[y$chr=="X"],Ne = my.Ne*2,t = c(0,my.gen))
  }
  covBase<-rpois(nrow(simTraj),Bcov)
  cntBase<-rbinom(nrow(simTraj),size = covBase,prob = simTraj[,1])
  covEvol<-rpois(nrow(simTraj),Ecov)
  cntEvol<-rbinom(nrow(simTraj),size=covEvol,prob=simTraj[,2])
  return(list(cntBase,covBase,cntEvol,covEvol))
}

#call.sims: call simulations
call.sims<-function(my.file,est_Ne,Bcov,Ecov){
  #obtain the generation 
  my.gen<-unlist(strsplit(my.file,"_"))
  my.gen<-my.gen[grep("F0",my.gen)]
  my.gen<-gsub("F0F","",my.gen)
  my.gen<-gsub("-ordered","",my.gen)
  my.gen<-as.numeric(my.gen)
  print(my.gen)
  my.data<-read.data(my.file)
  my.data<-calculate.average.ancestral(my.data)
  
  #simulations for autosomes
  print(Ecov)
  print(Bcov)
  sims.autosome<-replicate(10,calculate.counts(my.data,autosome = T,my.Ne = est_Ne[1],my.gen=my.gen,Bcov = Bcov,Ecov = Ecov),simplify = F)
  ListCnt<-sapply(sims.autosome,"[[",1,simplify = F)
  ListCnt<-append(ListCnt,sapply(sims.autosome,"[[",3,simplify = F))
  ListCov<-sapply(sims.autosome,"[[",2,simplify = F)
  ListCov<-append(ListCov,sapply(sims.autosome,"[[",4,simplify = F))
  snpNum<-length(ListCnt[[1]])
  print(length(ListCnt))
  #concatenate counts and coverages for all 10 replicates
  syncTable <- c()
  for (i in 1:length(ListCnt)) {
    syncTable <- cbind(syncTable, sapply(1:snpNum,function(x) {
      my.x<-ListCnt[[i]][x]
      my.y<-ListCov[[i]][x]
      convert.to.sync(x=my.x,y=my.y)
    }))
  }
  print(head(syncTable))
  base <- c(rep('A',time=snpNum))
  complete_syncTable <- cbind(as.character(my.data$chr[my.data$chr != 'X']),my.data$pos[my.data$chr != 'X'], base,syncTable)
  toSave<-paste("Dsim_Ne",est_Ne[1],"_autosome_",snpNum,"independentSNPs_10reps_simulated.sync",sep="")
  write.table(complete_syncTable,file = toSave,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  #simulations for X
  sims.X<-replicate(10,calculate.counts(my.data,autosome = F,my.Ne = est_Ne[2],my.gen=my.gen,Bcov = Bcov,Ecov = Ecov),simplify = F)
  ListCnt<-sapply(sims.X,"[[",1,simplify = F)
  ListCnt<-append(ListCnt,sapply(sims.X,"[[",3,simplify = F))
  ListCov<-sapply(sims.X,"[[",2,simplify = F)
  ListCov<-append(ListCov,sapply(sims.X,"[[",4,simplify = F))
  snpNum<-length(ListCov[[1]])
  #concatenate counts and coverages for all 10 replicates
  syncTable <- c()
  for (i in 1:length(ListCnt)) {
    syncTable <- cbind(syncTable, sapply(1:snpNum,function(x) {
      my.x<-ListCnt[[i]][x]
      my.y<-ListCov[[i]][x]
      convert.to.sync(x=my.x,y=my.y)
    }))
  }
  base <- c(rep('A',time=snpNum))
  complete_syncTable <- cbind(as.character(my.data$chr[my.data$chr == 'X']),my.data$pos[my.data$chr == 'X'], base,syncTable)
  toSave<-paste("Dsim_Ne",est_Ne[2],"_X_",snpNum,"independentSNPs_10reps_simulated.sync",sep="")
  write.table(complete_syncTable,file = toSave,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
}

########################
#ANALYSIS
########################

my.files<-list.files(pattern=glob2rx("*cov.txt"),path = "../data_emp",full.names = T)
print(my.files)

#form 02_Ne_estimation.R
my.ne<-readRDS("../data_emp/NeEstimates.rds")

#ancestral coverage based on Barghi et al 
neda_coverage<-read.csv("../data_emp/data_averageCoverage_PLoS.csv",header=T,sep=";",dec=",")
my_coverage<-ddply(neda_coverage,.(Generation),summarize,Coverage=mean(Average))
my_coverage$Coverage<-round(my_coverage$Coverage)

Bcov<-my_coverage$Coverage[my_coverage$Generation==0]

#F60
Ecov<-my_coverage$Coverage[my_coverage$Generation==60]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F60R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F60",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)

#F50
Ecov<-my_coverage$Coverage[my_coverage$Generation==50]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F50R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F50",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)

#F40
Ecov<-my_coverage$Coverage[my_coverage$Generation==40]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F40R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F40",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)

#F30
Ecov<-my_coverage$Coverage[my_coverage$Generation==30]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F30R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F30",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)

#F20
Ecov<-my_coverage$Coverage[my_coverage$Generation==20]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F20R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F20",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)

#F10
Ecov<-my_coverage$Coverage[my_coverage$Generation==10]
est_Ne<-my.ne[which(rownames(my.ne)%in%paste("F10R",c(1:10),sep="")),]
est_Ne<-as.numeric(round(colMeans(est_Ne)))
call.sims(my.file=my.files[grep("F10",my.files)],est_Ne=est_Ne,Bcov=Bcov,Ecov=Ecov)