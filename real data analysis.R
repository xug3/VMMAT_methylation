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


pheno=read.table("minfi/MESA_methy_all_no_dups_nomed_bmi.txt",header=T)
pheno[1,]
y.long=pheno$sbp
time=cbind(pheno$ID,pheno$age)
y.cov=pheno[,c(4:12,19,25)]
phe.model="Gaussian"
vmmat.est=vmmat_est(y.long,time,y.cov, phe.model=phe.model,VCcorrection=F)
vmmat.est$tau
vmmat.est$n.iter
glmm.est=glmm_est(y.long,time,y.cov,phe.model=phe.model)
glmm.est$tau
glmm.est$n.iter

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

  args = commandArgs(trailingOnly=TRUE)
  args= as.numeric(args[1])

  ncpg=dim(ADNI_methy)[1]
  nsep=ceiling(ncpg/100)

  if((args+(nsep-1)*100)<=ncpg){
    var.index=args+(0:(nsep-1))*100
  } else{
    var.index=args+(0:(nsep-2))*100
  }


  for(m in var.index){
    Gi <- c(unlist(ADNI_methy[m,]))
      VMMAT=vmmat_test(vmmat.est,Gi)
      GMMAT=glmm_test(glmm.est,Gi)[,2]
    pvalue=c(VMMAT,GMMAT)
    pvalue=as.matrix(t(pvalue),nrow=1)
    rownames(pvalue)=cpg.name[m]
    results=rbind(results,pvalue)
  }

colnames(results)=c("VMMAT","GMMAT")
write.table(results,paste0("zmethyall_",args,".txt"))
q("no")
