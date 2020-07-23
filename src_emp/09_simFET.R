########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: neutral simulations for an empirical p-value distribution of the FET tests, empirical data 

#original code: Neda Barghi "Genetic redundancy fuels polygenic adaptation in Drosophila", adapted by Anna Maria Langmüller 
rm(list=ls())
options(stringsAsFactors = FALSE)
my_wd<-"/Volumes/Temp1/early_response_Dsim_FL/src/submission/concordance_ER/data_emp/"
setwd(my_wd)
########################
#PACKAGES
########################
library(data.table)

########################
#FUNCTIONS
########################
#read.data: read tab separated data with header
read.data<-function(x,header=TRUE){
  initial<-read.table(x,header=header,sep="\t",nrows=100)
  classes<-sapply(initial,class)
  y<-read.table(x,header=header,sep="\t",colClasses = classes)
  return(y)
}


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

#calculate.counts: add noise onto allele counts, replicate specific
calculate.counts<-function(y,autosome,my.Ne,my.gen,my.rep,f.Bcov,f.Ecov){
  togrep<-paste("F0R",my.rep,sep="")
  my.subset<-c("chr","pos",togrep)
  temp.col<-y[,which(colnames(y)%in%my.subset)]
  if(autosome){
    print(paste("Launching Sim with Ne:",my.Ne,", autosome, Bcov",f.Bcov, "Ecov,",f.Ecov))
    simTraj<-fast.trajectory(p0=temp.col[temp.col$chr!="X",3],Ne=my.Ne*2,t=c(0,my.gen))
  } else {
    print(paste("Launching Sim with Ne:",my.Ne,",X, Bcov",f.Bcov, "Ecov,",f.Ecov))
    simTraj<-fast.trajectory(p0=temp.col[temp.col$chr=="X",3],Ne=my.Ne*2,t=c(0,my.gen))
  }
  covBase<-rpois(nrow(simTraj),f.Bcov)
  cntBase<-rbinom(nrow(simTraj),size = covBase,prob = simTraj[,1])
  covEvol<-rpois(nrow(simTraj),f.Ecov)
  cntEvol<-rbinom(nrow(simTraj),size=covEvol,prob=simTraj[,2])
  return(list(cntBase,covBase,cntEvol,covEvol))
}

#call.sims: call simulations
call.sims<-function(my.file,est_Ne,ff.Bcov,ff.Ecov,est_Ne_X,meanX=NULL,my.g){
  my.gen<-my.g
  my.data<-read.data(my.file)
  
  #simulations for autosomes
  sims.autosome<-list()
  sims.X<-list()
  #for each replicate
  for(i in c(1:10)){
    print(ff.Bcov[i])
    print(ff.Ecov[i])
    sims.autosome[[i]]<-calculate.counts(my.data,autosome = T,my.Ne = est_Ne[i],my.gen = my.gen,f.Bcov = ff.Bcov[i],f.Ecov = ff.Ecov[i],my.rep =i)
    print(str(sims.autosome))
    ListCnt<-list()
    ListCov<-list()
    ListCnt[[1]]<-sims.autosome[[i]][[1]]
    ListCnt[[2]]<-sims.autosome[[i]][[3]]
    ListCov[[1]]<-sims.autosome[[i]][[2]]
    ListCov[[2]]<-sims.autosome[[i]][[4]]
    snpNum<-length(ListCnt[[1]])
    print(snpNum)
    base <- c(rep('A',time=snpNum))
    syncTable<-c()
    for (j in 1:2) {
      syncTable <- cbind(syncTable, sapply(1:snpNum,function(x) {
        my.x<-ListCnt[[j]][x]
        my.y<-ListCov[[j]][x]
        convert.to.sync(x=my.x,y=my.y)
      }))
    }
    complete_syncTable <- cbind(as.character(my.data$chr[my.data$chr != 'X']),my.data$pos[my.data$chr != 'X'], base,syncTable)
    toSave<-paste("Dsim_Ne",est_Ne[i],"_autosome_",snpNum,"independentSNPs_rep",i,"F",my.gen,"_simulated.sync",sep="")
    write.table(complete_syncTable,file = toSave,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  
    if(any(my.data$chr=="X")){
    #simulations for X
	  if(is.null(meanX)){
    sims.X[[i]]<-calculate.counts(my.data,autosome = F,my.Ne = est_Ne_X[i],my.gen=my.gen,f.Bcov = ff.Bcov[i],f.Ecov = ff.Ecov[i],my.rep = i)
    } else {
        mean.Ne<-round(mean(est_Ne_X,na.rm = T))
        print(paste("Used mean Ne:",mean.Ne))
	    sims.X[[i]]<-calculate.counts(my.data,autosome = F,my.Ne = mean.Ne,my.gen = my.gen,f.Bcov = ff.Bcov[i],f.Ecov = ff.Ecov[i],my.rep = i)
    } 
    print(str(sims.X[[i]]))
    ListCnt<-list()
    ListCov<-list()
    ListCnt[[1]]<-sims.X[[i]][[1]]
    ListCnt[[2]]<-sims.X[[i]][[3]]
    ListCov[[1]]<-sims.X[[i]][[2]]
    ListCov[[2]]<-sims.X[[i]][[4]]
    snpNum<-length(ListCov[[1]])
    base <- c(rep('A',time=snpNum))
    syncTable <- c()
    for (j in 1:2) {
    syncTable <- cbind(syncTable, sapply(1:snpNum,function(x) {
      my.x<-ListCnt[[j]][x]
      my.y<-ListCov[[j]][x]
      convert.to.sync(x=my.x,y=my.y)
    }))
  }
  complete_syncTable <- cbind(as.character(my.data$chr[my.data$chr == 'X']),my.data$pos[my.data$chr == 'X'], base,syncTable)
  
  if(is.null(meanX)){
  toSave<-paste("Dsim_Ne",est_Ne_X[i],"_X_",snpNum,"independentSNPs_rep",i,"F",my.gen,"_simulated.sync",sep="")
  } else {
    toSave<-paste("Dsim_Ne",mean.Ne,"_MeanX_",snpNum,"independentSNPs_rep",i,"F",my.gen,"_simulated.sync",sep="")
  }
  write.table(complete_syncTable,file = toSave,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  } else {
      print("Skip!")
    }
}}


##################################
#ANALYSIS
##################################

my.coverages<-read.csv("data_averageCoverage_PLoS.csv",sep=";",dec=",")
my.nes<-readRDS("NeEstimates.rds")
my.files<-list.files(pattern="*cov.txt")
my.files<-rev(my.files)
print(my.files)
print(my.nes)
print(my.coverages)
my.generation<-seq(60,10,-10)
Bcov<-round(my.coverages$Average[which(my.coverages$Generation==0)])

#F60
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[1],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[1],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[1])])
print(Ecov)
call.sims(my.file = my.files[1],est_Ne = est_Ne_autosome ,ff.Bcov =Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[1])

#F50
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[2],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[2],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[2])])
print(Ecov)
call.sims(my.file = my.files[2],est_Ne = est_Ne_autosome ,ff.Bcov=Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[2])

#F40
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[3],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[3],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[3])])
print(Ecov)
call.sims(my.file = my.files[3],est_Ne = est_Ne_autosome ,ff.Bcov=Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[3])

#F30
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[4],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[4],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[4])])
print(Ecov)
call.sims(my.file = my.files[4],est_Ne = est_Ne_autosome ,ff.Bcov=Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[4])

#F20
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[5],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[5],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[5])])
print(Ecov)
call.sims(my.file = my.files[5],est_Ne = est_Ne_autosome ,ff.Bcov=Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[5])

#F10
est_Ne_autosome<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[6],"R",c(1:10),sep="")),1]
est_Ne_autosome<-round(est_Ne_autosome)
est_Ne_X<-my.nes[which(rownames(my.nes)%in%paste("F",my.generation[6],"R",c(1:10),sep="")),2]
Ecov<-round(my.coverages$Average[which(my.coverages$Generation==my.generation[6])])
print(Ecov)
call.sims(my.file = my.files[6],est_Ne = est_Ne_autosome ,ff.Bcov=Bcov,ff.Ecov=Ecov,est_Ne_X = est_Ne_X,meanX = T,my.g=my.generation[6])
