########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: haplotype file creation for forward simulations

########################
#PACKAGES
########################

rm(list=ls())
library(stringi)
library(stringr)
library(parallel)

########################
#FUNCTIONS
########################

#haplotype-file format for mimicrEE 2 :
#CHR POS REFERENCE ANCESTRAL/DERIVED HAPLOTYPE-INFORMATION 

derive.derived<-function(reference,haplotypes){
  temp<-table(unlist(strsplit(haplotypes,"")))
  temp<-temp[-which(names(temp)==" ")]
  if (any(names(temp)==as.character(reference))) temp<-temp[-which(names(temp)==as.character(reference))]
  output<-paste(reference,"/",names(temp[which.max(temp)]),sep="")
  return(output)
}


get.haplotypeInfo<-function(l.haplos,l.haplos.names,all.heterozygous=F,l.freqs=0,N=300,l.pick="Florida_155"){
  haplotypeInfo<-subset(l.haplos,select=c(chr,pos,M252))
  
  if (all.heterozygous & length(l.haplos.names)!=2)  {
    return(c("You can only get all heterozygous if you have only 2 haplotypes"))}
  
  if (all.heterozygous & length(l.haplos.names)==2) {
    print(paste("Creating heterozygous from ",l.haplos.names[1]," and ",l.haplos.names[2],sep=""))
    haplo.hetero<-paste(l.haplos[,grep(paste(l.haplos.names[1],"$",sep=""),colnames(l.haplos),perl=TRUE)],l.haplos[,grep(paste(l.haplos.names[2],"$",sep=""),colnames(l.haplos),perl=TRUE)],sep="")
    haplotypeInfo$haplotypes<-haplo.hetero
    
    for (i in c(2:N)){
      haplotypeInfo$haplotypes<-paste(haplotypeInfo$haplotypes,haplo.hetero,sep=" ")
    }
  }
  
  if (!all.heterozygous) {
    if (length(l.freqs)!=length(l.haplos.names)){return("You have to provide starting frequencies for all haplotypes!")} else {
      if(sum(l.freqs) > 1) return("The sum of your frequencies exceeds 1!")
      #print (paste("Creating homozygous of",l.haplos.names,"with frequencies",l.freqs,sep=" "))
      haplo.homo<-c()
      #get for main selected one a ceiling + create vector N.haplos that contains absolute number for the various haplotypes:
      N.haplos<-numeric(length(l.freqs))
      id.topick<-grep(paste(l.pick,"$",sep=""),l.haplos.names,perl = TRUE)
      N.haplos[id.topick]<-ceiling(N*l.freqs[id.topick])
      print(paste("Doing",N.haplos[id.topick],"Individuals for",l.haplos.names[id.topick],sep=" "))
      #for all others:
      for (f in c(1:length(l.freqs))[-id.topick]){
        rep.level<-floor(N*l.freqs[f])
        print(paste("Doing",rep.level,"Individuals for",l.haplos.names[f],sep=" "))
        N.haplos[f]<-rep.level}
      
      #if N.haplos is smaller than N- add random haplotypes to sum up:
      if(sum(N.haplos)<N) {
        to.add<-N-sum(N.haplos)
        set.seed(to.add)
        while(to.add>0){
          id.add<-sample(x = c(1:length(N.haplos)),size = 1)
          N.haplos[id.add]<-N.haplos[id.add]+1
          to.add<-to.add-1
        }
      }
      
      names(N.haplos)<-l.haplos.names
      print(N.haplos)
      print(sum(N.haplos))
      #do the actual haplotypes:
      haplo.homo<-c()
      for (f in c(1:length(l.haplos.names))){
        if(N.haplos[l.haplos.names[f]]!=0){
          temp.diploid<-paste(l.haplos[,l.haplos.names[f]],l.haplos[,l.haplos.names[f]],sep="")
          #print(paste("haplotype",f,"genotype",temp.diploid,"r:",N.haplos[l.haplos.names[f]],sep=" "))
          for (iter in c(1:N.haplos[l.haplos.names[f]])){
            haplo.homo<-paste(haplo.homo,temp.diploid,sep=" ")
          }}
      }
      #print(haplo.homo)
      haplotypeInfo$haplotypes<-haplo.homo 
      haplotypeInfo$ancestral.derived<-mapply(FUN = derive.derived,reference=haplotypeInfo$M252,haplotypes=haplotypeInfo$haplotypes)
      haplotypeInfo<-haplotypeInfo[c("chr","pos","M252","ancestral.derived","haplotypes")]
      #filter out monomorphic sites:
      haplotypeInfo<-subset(haplotypeInfo,nchar(ancestral.derived)==3)
      haplotypeInfo$haplotypes<-trimws(haplotypeInfo$haplotypes,which="both")
      return(haplotypeInfo)
    }
  }
}


is.polymorphic<-function(haplo.df,l.N=300){
  n.count<-rep(0,4)
  n.char<-c("A","T","C","G")
  for (n.char.iter in c(1:length(n.char))){
    n.count[n.char.iter]<-stri_count_fixed(haplo.df,n.char[n.char.iter])
  }
  if (any(n.count==2*l.N)| sum(n.count)!=2*l.N) {
    return (FALSE) } else {
      return(TRUE)
    }
}

#function for removing sites with InDels: 
is.indel<-function(haplo.df){
  single.haplos<-strsplit(as.character(haplo.df)," ")
  res<-any(unlist(lapply(single.haplos,nchar))!=2)
  return(res)
}

get.n<-function(x) as.numeric(gsub("-N","",str_extract(x,"-N[0-9]*")))

processFile = function(filepath,l.N=300) {
  con = file(filepath, "r")
  n<-get.n(filepath)
  print(paste("N:",n))
  l<-as.numeric(system(paste("cat",filepath,"| wc -l",sep=" "),intern = T))
  print(l)
  id.indel<-logical(length=l)
  id<-logical(length=l)
  counter<-1
  for (i in 1:l){
    linn=readLines(con,1)
    lin.f<-unlist(strsplit(linn,"\t"))
    #print(str(lin.f))
    id.indel[i]<-is.indel(lin.f[5])
    id[i]<-is.polymorphic(lin.f[5],l.N=n)
    if(i%%10000==0) {
      print(paste("Done:",round(counter*1000000/l,2),"%",sep=" "))
      counter<-counter+1
    }
  }
  close(con)
  print("Parsing done ...")
  sel<-c(id & !id.indel)
  print(table(sel))
  print("Id determination done...")
  dat<-read.table(filepath,header = F,sep = "\t",stringsAsFactors = F)
  print("Reading done ...")
  write.table(dat[sel,],file = gsub(".tab",".modified.tab",filepath),quote = F,row.names = F,col.names = F,sep = "\t")
  print("Writing done")
  return(c(l-length(which(sel))))
}

########################
#DATA & ANALYSIS
########################

#get the 189 available ancestral Florida haplotype files: 
base<-readRDS("../data_sim/all_biallelicSNPs_q50.rds")
print(head(base))
#after E-mail correspondence with Neda Barghi, 22/07/2019
id_Florida160<-which(colnames(base)=="other_17")
id_Florida184<-which(colnames(base)=="other_18")
colnames(base)[id_Florida160]<-"Florida_160"
colnames(base)[id_Florida184]<-"Florida_184"
rm(id_Florida160,id_Florida184)

id_Florida<-grep("Florida",colnames(base))
length(id_Florida)
#all 189 available
selected.haplos<-colnames(base)[id_Florida]

#generate MimicrEE2 file 
f<-rep(c(1/length(selected.haplos)),length(selected.haplos))
a<-get.haplotypeInfo(l.haplos = base,l.haplos.names = c(selected.haplos),all.heterozygous = F,l.freqs = f,l.pick = selected.haplos[1],N = 300)
write.table(x=a,file = "../data_sim/Dsim_FL_189Ancestral_wholeGenome-N300.tab",sep="\t",quote = F,row.names = F,col.names = F)
rm(a)
gc()
print("Done!")
my.files<-"../data_sim/Dsim_FL_189Ancestral_wholeGenome-N300.tab"
print(my.files)
processFile(filepath = my.files)
