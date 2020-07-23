########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Simulation without linkage, "winners curse"

rm(list=ls())
options(stringsAsFactors = F,scipen = 999)

library(poolSeq)
library(ggplot2)
library(ggpubr)

#Functions
#-------------------
compute_FDR <- function(emp_pval,sim_pval,thresh){ #Barghi, 2019 
  ecdf_sim <- ecdf(sim_pval)
  ecdf_emp <- ecdf(emp_pval)
  fdr <- (1-ecdf_sim(emp_pval)) / (1-ecdf_emp(emp_pval))
  fdr[emp_pval > max(sim_pval)] <- 0
  fdr[sim_pval < min(emp_pval)] <- 1
  candidates <- emp_pval[fdr <= thresh]
  FDR <- min(candidates, na.rm = TRUE)
  return(FDR)
}

#x ... matrix with simulated trajectories (ncol = 7, nrow = nr * ns)
#g ... generations to look at (length = 2)
#nr ... number of replicates 
#ns ... number of sites 
compute_CMH<-function(x,g,nr,ns,ncov=100){
  A0<-x[,c(g[1]/10+1)]
  A0<-matrix(A0,nrow = nr,ncol = ns)
  A0<-round(A0*ncov)
  a0<-ncov-A0
  At<-x[,c(g[2]/10+1)]
  At<-matrix(At,nrow = nr,ncol = ns)
  At<-round(At*ncov)
  at<-ncov-At
  y<-cmh.test(A0=A0,a0=a0,At=At,at=at,log = TRUE)
  return(y)
}

#x .... simulation matrix
#g ... generations to look at 
#nr ... number of replicates 
#ns ... number of snps simulated
#thres ... min frequency change
compute_RISING<-function(x,g,nr,ns,thres){
  A1<-x[,c(g[1]/10+1)]
  A2<-x[,c(g[2]/10+1)]
  A1<-matrix(A1,nrow = nr,ncol = ns)
  A2<-matrix(A2,nrow = nr,ncol = ns)
  dA<-A2-A1
  reps<-apply(dA,2,function(x) sum(x>=thres))
  return(reps)
}
#-------------------
#Variables:
s<-0.05 #median s Barghi et al. (2019)
saf<-0.10 #median starting allele frequency Barghi et al. (2019)
ne<-300
h<-0.5 
gens<-seq(0,60,10)
nreps<-10
nsnp<-100000
my_cov<-100
p.value<-c()
for(nthres in c(0.2,0.15,0.1,0.05)){
  sim<-wf.traj(p0 = rep(saf,nreps*nsnp),Ne = ne,t = gens,s = s,h = h,haploid = FALSE,approximate = FALSE)
  sim_neutral<-wf.traj(p0 = rep(saf,nreps*nsnp),Ne = ne,t = gens,s=0,h=h,haploid = FALSE,approximate = FALSE)
  cmh_sim<-compute_CMH(sim,c(0,20),nreps,nsnp,my_cov)
  cmh_neutral<-compute_CMH(sim_neutral,c(0,20),nreps,nsnp,my_cov)

  print(summary(cmh_neutral))
  print(summary(cmh_sim))

  fdr<-compute_FDR(emp_pval = cmh_sim,sim_pval = cmh_neutral,thresh = 0.05)
  print(fdr)
  print(table(cmh_sim>=fdr))

  rising<-compute_RISING(sim,gens,nreps,nsnp,nthres)
  n_df<-data.frame(nrep=rising,detected=cmh_sim>=fdr)
  table(n_df$detected)
  p.value<-c(p.value,wilcox.test(n_df$nrep~n_df$detected)$p.value)
  print(p.value)
}


p.value<-p.adjust(p.value,method = "BH")
cols <- c("FALSE" = "darkorchid4", "TRUE" = "cornflowerblue")
n_df$detected<-factor(n_df$detected,levels = c("TRUE","FALSE"),ordered = T)

g<-ggboxplot(data=n_df,x="detected",y="nrep",fill="detected")+scale_fill_manual(name="generation 20",values=cols,breaks=c("TRUE","FALSE"),labels=c("candidate", "no-candidate"))+scale_x_discrete("generation 20",labels=c("candidate","no candidate"))
g<-g+scale_y_continuous("number of rising replicates in generation 20",breaks=c(0,2,4,6,8,10),labels=c("0","2","4","6","8","10"))
g<-g+coord_cartesian(ylim=c(0,14))
p_af<-data.frame(group1=c("TRUE"),group2=c("FALSE"),p.adj=round(p.value[4],4),y.position=c(12),x.position=c(1.5),stringsAsFactors = F)
g<-g+stat_pvalue_manual(p_af,label = "p.adj")

g
print(round(p.value,3))
