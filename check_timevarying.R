library(data.table)
library(MASS)
library(nleqslv)
library(mvtnorm,lib)
library(expm)
library(CompQuadForm)
library(readxl)


source("vmmat_est.R")
source("vmmat_util.R")
source("vmmat_pvalue.R")


pheno=read.table("minfi/MESA_methy_all_no_dups_nomed.txt",header=T)
pheno[1:10,]
y.long=pheno$sbp

time=cbind(pheno$ID,pheno$age)
y.cov=pheno[,c(4:12,19)]


results=NULL

ADNI_methy=fread(paste0("minfi0.01/raw_m.txt"),header=F)
ADNI_methy_header=fread("minfi0.01/header_file.txt",header=F)
ADNI_methy_header=unlist(ADNI_methy_header)


ADNI_methy=data.frame(ADNI_methy)
rownames(ADNI_methy)=ADNI_methy[,1]
ADNI_methy=ADNI_methy[,-1]
ADNI_methy=ADNI_methy[,match(pheno$Sample_ID,ADNI_methy_header)]

dim( ADNI_methy)
cpg.name=rownames(ADNI_methy)



table(time[,2])
effect_all=matrix(0,5,7)
rownames(effect_all)=c("cg08196267","cg08023851","cg27243121","cg03475420","cg05550588")

var.index=c(which(cpg.name=="cg08196267"),which(cpg.name=="cg08023851"),which(cpg.name=="cg27243121"),which(cpg.name=="cg03475420"),which(cpg.name=="cg05550588"))


samplesize=rep(0,7)
results=NULL
for(k in 1:5){
  m=var.index[k]
  Gi <- c(unlist(ADNI_methy[m,]))
  data_all=cbind(y.long,time,y.cov,Gi)
  table(data_all[,3])
  data_alli=data_all[data_all[,3]<55,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,1]=modeli$coefficients[12]
  samplesize[1]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=55&data_all[,3]<60,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,2]=modeli$coefficients[12]
  samplesize[2]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=60&data_all[,3]<65,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,3]=modeli$coefficients[12]
  samplesize[3]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=65&data_all[,3]<70,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,4]=modeli$coefficients[12]
  samplesize[4]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=70&data_all[,3]<75,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,5]=modeli$coefficients[12]
  samplesize[5]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=75&data_all[,3]<80,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,6]=modeli$coefficients[12]
  samplesize[6]=dim(data_alli)[1]
  
  data_alli=data_all[data_all[,3]>=80,-c(2,3)]
  dim(data_alli)
  modeli=lm(y.long~.,data=data_alli)
  effect_all[k,7]=modeli$coefficients[12]
  samplesize[7]=dim(data_alli)[1]
  
}


plot(effect_all[1,], xaxt = "n",type="o",ylim=c(min(effect_all),max(effect_all)),
     ylab="Regression coefficients",xlab="",main="Methlytation effects on SBP in different age group")
lines(effect_all[2,],type="o",col=2)
lines(effect_all[3,],type="o",col=3)
lines(effect_all[4,],type="o",col=4)
lines(effect_all[5,],type="o",col=5)
axis(1, at=1:7, labels=c("44-55","55-60","60-65","65-70","70-75","75-80","80-93"))
axis(1, at = 1:7, labels = samplesize, line = 1.2, tick = FALSE)
mtext("N_obs",     side = 1, line = 2.1, at = 0.4, adj = 0)

legend("topright",legend = rownames(effect_all),col    = 1:5, lty    = 1,                     # solid line
  lwd    = 2,pch    = 1, bty    = "n" )

dev.off()
